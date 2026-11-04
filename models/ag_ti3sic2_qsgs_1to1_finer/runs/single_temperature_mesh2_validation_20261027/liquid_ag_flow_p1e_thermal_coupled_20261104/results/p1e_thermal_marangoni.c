/*
 * P1e numerical integration test: CLSVOF + integral.h with an advected,
 * diffused dimensionless temperature tracer driving sigma(T).
 * No Ag properties, material geometry, or production thermal solver.
 */
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "diffusion.h"
#include "integral.h"

int cells_per_radius = 16, LEVEL = 8;
int use_amr = 1;
double PeT = 0.01, gradT = 1., t_end, Udrop, initial_mass, xcm0, tprev;
double end_t0 = 3.;
const double R = 1. [1], Mu = 1., Rho = 1. [0];
const double Re = 0.066, Ca = 0.066;
const double GammaT = Re*sq(Mu)/(Rho*sq(R));
const double Gamma0 = GammaT*R/Ca;
const double t0 = Mu/(GammaT);
const double Uref = GammaT*R/Mu;
double thermal_diffusivity;
scalar Tfield[], sigmaf[];
face vector thermal_D[];

Tfield[left] = dirichlet(-gradT*L0/2.);
Tfield[right] = dirichlet(gradT*L0/2.);
Tfield[bottom] = neumann(0.);
Tfield[top] = neumann(0.);

int main (int argc, char * argv[])
{
  if (argc < 4 || argc > 6) {
    fprintf (stderr, "usage: %s cells_per_radius PeT gradient [end_time/t0] [use_amr: 0|1]\n", argv[0]);
    return 2;
  }
  cells_per_radius = atoi(argv[1]);
  PeT = atof(argv[2]);
  gradT = atof(argv[3]);
  if (argc >= 5) end_t0 = atof(argv[4]);
  if (argc >= 6) use_amr = atoi(argv[5]);
  if ((cells_per_radius != 8 && cells_per_radius != 16 && cells_per_radius != 32) ||
      PeT <= 0. || (gradT != 0. && gradT != 1.) || end_t0 <= 0. ||
      (use_amr != 0 && use_amr != 1))
    return 2;
  thermal_diffusivity = Uref*R/PeT;
  LEVEL = (int)(log(16.*cells_per_radius)/log(2.) + 0.5);
  N = 16*cells_per_radius;
  size (16.*R);
  origin (-L0/2.);
  rho1 = rho2 = Rho;
  mu1 = mu2 = Mu;
  d.sigmaf = sigmaf;
  tracers = list_concat (tracers, {Tfield});
  TOLERANCE = 1e-4 [*];
  Udrop = -2./((2. + 3.*mu2/mu1)*(2. + 1.))*GammaT*R*gradT/mu1;
  t_end = (gradT == 0. ? 0.5 : end_t0)*t0;
  run();
}

/* Init the passive temperature with the steady linear conduction solution. */
event init (t = 0)
{
  foreach() {
    d[] = sqrt(sq(x) + sq(y)) - R;
    Tfield[] = gradT*x;
    sigmaf[] = Gamma0 + GammaT*Tfield[];
  }
  boundary ({d,Tfield,sigmaf});
  foreach_face() {
    thermal_D.x[] = thermal_diffusivity*fm.x[];
  }
}

/* Centered solver transports Tfield through the tracer list; this implicit
   diffusion event restores a prescribed axial gradient at the Dirichlet ends.
   sigmaf is refreshed after thermal transport and before the next force update. */
event tracer_diffusion (i++)
{
  diffusion (Tfield, dt, thermal_D);
  foreach()
    sigmaf[] = Gamma0 + GammaT*Tfield[];
  boundary ({Tfield,sigmaf});
}

double max_u = 0., interface_T_rms = 0., measured_speed = 0.;
event logfile (i += 10)
{
  double xb=0., vb=0., mass=0., umax=0.;
  double wsum=0., terr2=0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:mass)
          reduction(max:umax) reduction(+:wsum) reduction(+:terr2)) {
    double v=(1.-f[])*dv();
    xb += x*v;
    vb += u.x[]*v;
    mass += v;
    umax=max(umax,norm(u));
    double gx=(f[1]-f[-1])/(2.*Delta);
    double gy=(f[0,1]-f[0,-1])/(2.*Delta);
    double w=sqrt(sq(gx)+sq(gy))*dv();
    wsum += w;
    terr2 += w*sq(Tfield[]-gradT*x);
  }
  double xc=mass>0. ? xb/mass : 0.;
  measured_speed=t>tprev ? (xc-xcm0)/(t-tprev) : 0.;
  max_u=umax;
  interface_T_rms=wsum>0. ? sqrt(terr2/wsum) : 0.;
  if (t==0.) {
    initial_mass=mass;
    xcm0=xc;
    tprev=0.;
    fprintf(stdout,"t_over_t0 mass_drift_frac xcm mean_u_over_Uref interval_speed_over_Uref umax interface_T_rms dt cells mgp_i mgp_resa PeT gradT\n");
  }
  fprintf(stdout,"%.9g %.12g %.12g %.12g %.12g %.12g %.12g %.9g %ld %d %.12g %.6g %.6g\n",
          t/t0,initial_mass>0.?(mass-initial_mass)/initial_mass:0.,xc,
          vb/mass/Uref,measured_speed/Uref,umax,interface_T_rms,dt,grid->tn,
          mgp.i,mgp.resa,PeT,gradT);
  fflush(stdout);
  xcm0=xc;
  tprev=t;
}

event report (t = t_end)
{
  double mass=statsf(f).sum;
  fprintf(stderr,"REPORT cells_per_radius=%d PeT=%.8g gradT=%.8g measured_speed_over_Uref=%.12g theory_over_Uref=%.12g rel_error_pct=%.12g max_u=%.12g interface_T_rms=%.12g mass=%.12g t_over_t0=%.12g cells=%ld mgp_i=%d mgp_resa=%.12g\n",
          cells_per_radius,PeT,gradT,measured_speed/Uref,Udrop/Uref,
          gradT ? 100.*fabs(measured_speed/Udrop-1.) : 0.,max_u,interface_T_rms,
          mass,t_end/t0,grid->tn,mgp.i,mgp.resa);
}

event stop (t = t_end)
{
  return 1;
}

#if TREE
event adapt (i++)
{
  if (use_amr)
    adapt_wavelet ({f,u,Tfield}, (double[]){1e-2,1e-5,1e-5,1e-5}, LEVEL);
}
#endif
