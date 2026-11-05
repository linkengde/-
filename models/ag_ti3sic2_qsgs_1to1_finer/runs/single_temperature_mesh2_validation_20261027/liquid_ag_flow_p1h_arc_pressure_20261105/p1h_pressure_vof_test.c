#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "iforce.h"
#include "curvature.h"
#include "p1h_arc_pressure.h"

#define TEST_N 32

static scalar kappa_field[], phi_pressure[];
static scalar pressure_effective_field[];
static P1HArcPressure pressure_cfg;
static double column_height[TEST_N*TEST_N];
static const double sigma_test = 1.;
static const double z_initial = 0.62;
static double test_dt = 1.e-4;
static double test_end = 0.2;
static double test_pressure = 1.67;
static double volume_initial = 0.;
static int local_test = 0;

u.n[back] = dirichlet(0.);
u.t[back] = neumann(0.);
u.r[back] = neumann(0.);
u.n[front] = dirichlet(0.);
u.t[front] = neumann(0.);
u.r[front] = neumann(0.);
p[front] = neumann(0.);

int main (int argc, char ** argv)
{
  if (argc != 2 || (strcmp(argv[1],"uniform") && strcmp(argv[1],"gaussian"))) {
    fprintf(stderr,"usage: %s uniform|gaussian\n",argv[0]);
    return 2;
  }
  local_test = !strcmp(argv[1],"gaussian");
  const char * value = getenv("P1H_TEST_DT");
  if (value) test_dt = strtod(value,NULL);
  value = getenv("P1H_TEST_END");
  if (value) test_end = strtod(value,NULL);
  value = getenv("P1H_TEST_PRESSURE");
  if (value) test_pressure = strtod(value,NULL);
  if (!(test_dt > 0.) || !(test_end > 0.) || !(test_pressure > 0.)) {
    fprintf(stderr,"invalid P1H test environment setting\n");
    return 2;
  }
  size(1.);
  origin(0.,0.,0.);
  N = TEST_N;
  periodic(left);
  periodic(bottom);
  rho1 = rho2 = 1.;
  mu1 = mu2 = 0.05;
  CFL = 0.2;
  DT = test_dt;
  TOLERANCE = 1.e-8;
  pressure_cfg = (P1HArcPressure){
    /* The Gaussian test uses Pi_p = p R / sigma ~= 0.30; this is a
     * dimensionless numerical benchmark, not an Ag arc-pressure estimate. */
    /* Pressure here is in the nondimensional benchmark's solver units, not Pa. */
    .enabled = 1, .peak_Pa = local_test ? test_pressure : 0.2,
    .center_x_m = 0.5, .center_y_m = 0.5, .radius_m = 0.18,
    .start_s = 0., .ramp_s = 0.002, .end_s = local_test ? test_end : fmin(test_end,0.04),
    .period_x_m = 1., .period_y_m = 1.,
    .profile = local_test ? P1H_PRESSURE_PERIODIC_GAUSSIAN : P1H_PRESSURE_UNIFORM
  };
  run();
}

event init (t = 0)
{
  vertex scalar plane[];
  foreach_vertex()
    plane[] = z_initial - z;
  fractions(plane,f);
  boundary({f});
  double initial_volume = 0.;
  foreach(reduction(+:initial_volume))
    initial_volume += f[]*dv();
  volume_initial = initial_volume;
}

event acceleration (i++)
{
  curvature(f,kappa_field);
  phi_pressure = new scalar;
  foreach() {
    double p_eff = 0.;
    double grad_z = 0.;
    if (f[] > 1.e-8 && f[] < 1. - 1.e-8) {
      double gx=(f[1]-f[-1])/(2.*Delta);
      double gy=(f[0,1]-f[0,-1])/(2.*Delta);
      grad_z=(f[0,0,1]-f[0,0,-1])/(2.*Delta);
      double gn=sqrt(gx*gx+gy*gy+grad_z*grad_z);
      if (gn > 1.e-12) {
        double nz=grad_z/gn;
        p_eff=p1h_arc_pressure_effective(&pressure_cfg,x,y,t,1.,1.,nz);
        phi_pressure[]=sigma_test*kappa_field[]+p_eff;
      }
      else
        phi_pressure[]=sigma_test*kappa_field[];
    }
    else {
      phi_pressure[]=nodata;
      p_eff=0.;
    }
    pressure_effective_field[]=p_eff;
  }
  boundary({phi_pressure,pressure_effective_field});
  double pressure_force_z = 0.;
  foreach_face(z, reduction(+:pressure_force_z)) {
    int current_interface=f[]>1.e-8 && f[]<1.-1.e-8;
    int lower_interface=f[0,0,-1]>1.e-8 && f[0,0,-1]<1.-1.e-8;
    if (f[] != f[0,0,-1] && (current_interface || lower_interface)) {
      double p_face=current_interface && lower_interface ?
        0.5*(pressure_effective_field[]+pressure_effective_field[0,0,-1]) :
        current_interface ? pressure_effective_field[] : pressure_effective_field[0,0,-1];
      pressure_force_z += p_face*(f[]-f[0,0,-1])/Delta*dv();
    }
  }
  f.phi=phi_pressure;
  if (i == 0)
    fprintf(stderr,"P1H_VOF_TEST mode=%s peak_test_pressure_units=%.9g dt_max=%.9g end_time=%.9g initial_volume=%.12g\n",
            local_test ? "periodic_gaussian" : "uniform",
            pressure_cfg.peak_Pa,test_dt,test_end,volume_initial);
  if (i%10 == 0)
    fprintf(stderr,"P1H_VOF_TEST_DIAG t=%.9g force_z=%.12g\n",t,pressure_force_z);
}

event log_state (i++)
{
  double volume=0., max_speed=0.;
  memset(column_height,0,sizeof(column_height));
  foreach(reduction(+:volume) reduction(max:max_speed)) {
    volume += f[]*dv();
    max_speed=max(max_speed,sqrt(sq(u.x[])+sq(u.y[])+sq(u.z[])));
  }
  foreach(serial) {
    int ix=clamp((int)(x*TEST_N),0,TEST_N-1);
    int iy=clamp((int)(y*TEST_N),0,TEST_N-1);
    column_height[ix*TEST_N+iy] += f[]*Delta;
  }
  double center_height=0., flank_height=0., center_count=0., flank_count=0.;
  for (int ix=0; ix<TEST_N; ix++)
    for (int iy=0; iy<TEST_N; iy++) {
      double xc=(ix+0.5)/TEST_N, yc=(iy+0.5)/TEST_N;
      double r2=sq(xc-0.5)+sq(yc-0.5);
      if (r2<sq(0.12)) { center_height+=column_height[ix*TEST_N+iy]; center_count++; }
      if (r2>sq(0.28)) { flank_height+=column_height[ix*TEST_N+iy]; flank_count++; }
    }
  if (center_count>0.) center_height/=center_count;
  if (flank_count>0.) flank_height/=flank_count;
  if (i%10 == 0)
    fprintf(stderr,"P1H_VOF_STATE t=%.9g center_h=%.12g flank_h=%.12g volume_rel=%.9g max_u=%.9g\n",
            t,center_height,flank_height,(volume-volume_initial)/volume_initial,max_speed);
}

event finish (t = test_end)
{
  double volume=0.;
  memset(column_height,0,sizeof(column_height));
  foreach(reduction(+:volume)) {
    volume += f[]*dv();
  }
  foreach(serial) {
    int ix=clamp((int)(x*TEST_N),0,TEST_N-1);
    int iy=clamp((int)(y*TEST_N),0,TEST_N-1);
    column_height[ix*TEST_N+iy] += f[]*Delta;
  }
  double center_height=0., flank_height=0., center_count=0., flank_count=0.;
  for (int ix=0; ix<TEST_N; ix++)
    for (int iy=0; iy<TEST_N; iy++) {
      double xc=(ix+0.5)/TEST_N, yc=(iy+0.5)/TEST_N;
      double r2=sq(xc-0.5)+sq(yc-0.5);
      if (r2<sq(0.12)) { center_height+=column_height[ix*TEST_N+iy]; center_count++; }
      if (r2>sq(0.28)) { flank_height+=column_height[ix*TEST_N+iy]; flank_count++; }
    }
  if (center_count>0.) center_height/=center_count;
  if (flank_count>0.) flank_height/=flank_count;
  double drift=fabs((volume-volume_initial)/volume_initial);
  double displacement=local_test ? flank_height-center_height : fabs(center_height-z_initial);
  double dx=1./TEST_N;
  printf("P1H_VOF_RESULT mode=%s center_height=%.12g flank_height=%.12g displacement=%.12g volume_drift=%.12g threshold=%.12g\n",
         local_test ? "gaussian" : "uniform",center_height,flank_height,
         displacement,drift,0.5*dx);
  if (drift > 1.e-5) {
    fprintf(stderr,"FAIL mass drift %.9g\n",drift);
    exit(1);
  }
  if (local_test && !(center_height < flank_height - 1.e-4)) {
    fprintf(stderr,"FAIL positive local pressure did not depress the center surface\n");
    exit(1);
  }
  if (!local_test && displacement > 0.5*dx) {
    fprintf(stderr,"FAIL uniform planar pressure changed the flat surface beyond half a cell\n");
    exit(1);
  }
  printf("PASS\n");
  exit(0);
}
