/*
 * P1F numerical sensitivity driver.  P1F_MODE=1 varies only DT with PeT fixed
 * at 0.01; P1F_MODE=2 varies only PeT with DT fixed at DT_BASE.
 * This is a dimensionless benchmark, not a material Ag model.
 */
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "diffusion.h"
#include "integral.h"

#define P1F_TIMESTEP 1
#define P1F_DIFFUSIVITY 2
#define P1F_STATIC_BALANCE 3
#define P1F_STATIC_DIAGNOSTIC 4

int cells_per_radius = 16;
int LEVEL = 8;
const double DT_BASE = 0.00506658613; /* largest logged P1E R16 baseline step */
const double PeT_FIXED = 0.01;
double gradT = 1.;
const double R = 1. [1], Mu = 1., Rho = 1. [0];
const double Re = 0.066, Ca = 0.066;
const double GammaT = Re*sq(Mu)/(Rho*sq(R));
const double Gamma0 = 1.;
const double t0 = Mu/(GammaT);
const double Uref = GammaT*R/Mu;
double end_t0 = 0.5;
const double fit_window_start_t0 = 0.30;

double PeT = PeT_FIXED, dt_factor = 1., thermal_diffusivity, t_end;
double tolerance_factor = 1.;
int static_pressure_balance = 0, static_pressure_sign = 1;
int use_amr = 1, amr_mode = 1;
scalar Tfield[], sigmaf[];
face vector thermal_D[];
FILE * history_file = NULL;
double initial_mass = 0., max_abs_mass_drift = 0.;
double last_max_velocity = 0.;

Tfield[left] = dirichlet(-gradT*L0/2.);
Tfield[right] = dirichlet(gradT*L0/2.);
Tfield[bottom] = neumann(0.);
Tfield[top] = neumann(0.);

int main (int argc, char * argv[])
{
#if P1F_MODE == P1F_TIMESTEP
  if (argc < 2 || argc > 3) {
    fprintf (stderr, "usage: %s dt_factor (1, 0.5, or 0.25) [gradient 0|1]\n", argv[0]);
    return 2;
  }
  dt_factor = atof(argv[1]);
  if (argc == 3) gradT = atof(argv[2]);
  if (dt_factor != 1. && dt_factor != 0.5 && dt_factor != 0.25)
    return 2;
  if (gradT != 0. && gradT != 1.)
    return 2;
  PeT = PeT_FIXED;
#elif P1F_MODE == P1F_DIFFUSIVITY
  if (argc < 2 || argc > 3) {
    fprintf (stderr, "usage: %s PeT (0.005, 0.01, 0.02, or 0.05) [gradient 0|1]\n", argv[0]);
    return 2;
  }
  PeT = atof(argv[1]);
  if (argc == 3) gradT = atof(argv[2]);
  if (PeT != 0.005 && PeT != 0.01 && PeT != 0.02 && PeT != 0.05)
    return 2;
  if (gradT != 0. && gradT != 1.)
    return 2;
  dt_factor = 1.;
#elif P1F_MODE == P1F_STATIC_BALANCE
  if (argc > 2) {
    fprintf (stderr, "usage: %s [pressure_sign: -1|1]\n", argv[0]);
    return 2;
  }
  if (argc == 2) static_pressure_sign = atoi(argv[1]);
  if (static_pressure_sign != -1 && static_pressure_sign != 1)
    return 2;
  PeT = PeT_FIXED;
  dt_factor = 0.25;
  gradT = 0.;
  static_pressure_balance = 1;
#elif P1F_MODE == P1F_STATIC_DIAGNOSTIC
  if (argc < 7 || argc > 8) {
    fprintf (stderr, "usage: %s cells_per_radius(16|32) amr_mode(0=uniform,1=dynamic,2=frozen,3=interface-only criteria) dt_factor pressure_sign(-1|0|1) end_t0 tolerance_factor(1|0.1|0.01) [gradient(0|1)]\n", argv[0]);
    return 2;
  }
  cells_per_radius = atoi(argv[1]);
  amr_mode = atoi(argv[2]);
  use_amr = amr_mode != 0;
  dt_factor = atof(argv[3]);
  static_pressure_sign = atoi(argv[4]);
  end_t0 = atof(argv[5]);
  tolerance_factor = atof(argv[6]);
  gradT = argc == 8 ? atof(argv[7]) : 0.;
  if ((cells_per_radius != 16 && cells_per_radius != 32) ||
      (amr_mode < 0 || amr_mode > 3) ||
      !(dt_factor == 0.25 || dt_factor == 0.125) ||
      (static_pressure_sign != -1 && static_pressure_sign != 0 && static_pressure_sign != 1) ||
      !(tolerance_factor == 1. || tolerance_factor == 0.1 || tolerance_factor == 0.01) ||
      (gradT != 0. && gradT != 1.) ||
      !(end_t0 >= 0.05 && end_t0 <= 0.5))
    return 2;
  PeT = PeT_FIXED;
  static_pressure_balance = static_pressure_sign != 0;
#endif
  DT = DT_BASE*dt_factor;
  thermal_diffusivity = Uref*R/PeT;
  t_end = end_t0*t0;
  N = 16*cells_per_radius;
  LEVEL = 0;
  for (int n = N; n > 1; n /= 2)
    LEVEL++;
  size (16.*R);
  origin (-L0/2.);
  rho1 = rho2 = Rho;
  mu1 = mu2 = Mu;
  d.sigmaf = sigmaf;
  /* The default tracer list is a static compound literal; copy before append. */
  tracers = list_concat (tracers, {Tfield});
  TOLERANCE = 1e-4 [*];
  TOLERANCE *= tolerance_factor;
  run();
}

event init (t = 0)
{
  foreach() {
    d[] = sqrt(sq(x) + sq(y)) - R;
    Tfield[] = gradT*x;
    sigmaf[] = Gamma0 + GammaT*Tfield[];
    if (static_pressure_balance)
      p[] = d[] < 0. ? static_pressure_sign*2.*Gamma0/R : 0.;
  }
  boundary ({d,Tfield,sigmaf,p});
  foreach_face()
    thermal_D.x[] = thermal_diffusivity*fm.x[];
  history_file = fopen("xcm_history.csv", "w");
  if (!history_file) {
    fprintf (stderr, "could not open xcm_history.csv for writing\n");
    exit (2);
  }
  fprintf (history_file,
           "t_over_t0,xcm,mass_drift_frac,max_velocity,interface_T_mean,"
           "interface_T_rms_departure,interface_gradT_tangent_rms,"
           "sigma_interface_min,sigma_interface_max,dt,cells,mgp_resa,"
           "liquid_u_rms,liquid_u_max,p_center,p_outer,pressure_jump,"
           "meridional_curvature_mean,meridional_curvature_rms,"
           "meridional_curvature_min,meridional_curvature_max,mgp_i,mgp_resb,"
           "mgp_resa_over_resb,mgpf_i,mgpf_resb,mgpf_resa,mgpf_resa_over_resb\n");
}

/* Advance temperature with the same velocity as the liquid and then refresh
   sigma(T) before the next interfacial force evaluation. */
event tracer_diffusion (i++)
{
  diffusion (Tfield, dt, thermal_D);
  foreach()
    sigmaf[] = Gamma0 + GammaT*Tfield[];
  boundary ({Tfield,sigmaf});
}

/* One record per physical timestep makes the velocity fit use the same time
   interval even when DT changes.  The LSQ fit is performed by the analyzer. */
event diagnostics (i++)
{
  double xb = 0., mass = 0., umax = 0., liquid_umax = 0., liquid_u2 = 0.;
  double wsum = 0., tsum = 0., tdev2 = 0., gt2 = 0.;
  double kap1 = 0., kap2 = 0., kapmin = HUGE, kapmax = -HUGE;
  double sigma_min = HUGE, sigma_max = -HUGE;
  foreach (reduction(+:xb) reduction(+:mass) reduction(max:umax)
           reduction(max:liquid_umax) reduction(+:liquid_u2)
           reduction(+:wsum) reduction(+:tsum) reduction(+:tdev2)
           reduction(+:gt2) reduction(+:kap1) reduction(+:kap2)
           reduction(min:kapmin) reduction(max:kapmax)
           reduction(min:sigma_min) reduction(max:sigma_max)) {
    double v = (1. - f[])*dv();
    mass += v;
    xb += x*v;
    umax = max(umax, norm(u));
    liquid_umax = max(liquid_umax, (1. - f[]) > 1e-6 ? norm(u) : 0.);
    liquid_u2 += v*(sq(u.x[]) + sq(u.y[]));
    double gx = (f[1] - f[-1])/(2.*Delta);
    double gy = (f[0,1] - f[0,-1])/(2.*Delta);
    double w = sqrt(sq(gx) + sq(gy))*dv();
    if (w > 0.) {
      double Tx = (Tfield[1] - Tfield[-1])/(2.*Delta);
      double Ty = (Tfield[0,1] - Tfield[0,-1])/(2.*Delta);
      double dx = (d[1] - d[-1])/(2.*Delta);
      double dy = (d[0,1] - d[0,-1])/(2.*Delta);
      double dn = sqrt(sq(dx) + sq(dy));
      double tangent_grad = 0.;
      if (dn > 0.) {
        double normal_grad = (Tx*dx + Ty*dy)/dn;
        tangent_grad = sqrt(max(0., sq(Tx) + sq(Ty) - sq(normal_grad)));
      }
      wsum += w;
      tsum += Tfield[]*w;
      tdev2 += sq(Tfield[] - gradT*x)*w;
      gt2 += sq(tangent_grad)*w;
    }
    if (fabs(d[]) < 2.*Delta) {
      sigma_min = min(sigma_min, sigmaf[]);
      sigma_max = max(sigma_max, sigmaf[]);
      double kap = distance_curvature (point, d);
      kap1 += kap*w;
      kap2 += sq(kap)*w;
      kapmin = min(kapmin, kap);
      kapmax = max(kapmax, kap);
    }
  }
  double xcm = mass > 0. ? xb/mass : 0.;
  last_max_velocity = umax;
  if (initial_mass == 0.)
    initial_mass = mass;
  double drift = initial_mass > 0. ? (mass - initial_mass)/initial_mass : 0.;
  max_abs_mass_drift = max(max_abs_mass_drift, fabs(drift));
  double interface_T_mean = wsum > 0. ? tsum/wsum : 0.;
  double interface_T_rms = wsum > 0. ? sqrt(tdev2/wsum) : 0.;
  double interface_gradT_tangent_rms = wsum > 0. ? sqrt(gt2/wsum) : 0.;
  double liquid_urms = mass > 0. ? sqrt(liquid_u2/mass) : 0.;
  double kapmean = wsum > 0. ? kap1/wsum : 0.;
  double kaprms = wsum > 0. ? sqrt(kap2/wsum) : 0.;
  if (sigma_min == HUGE) sigma_min = sigma_max = 0.;
  double pcenter = interpolate (p, 0., 0.);
  double pout = interpolate (p, 0., 1.5*R);
  double mgp_ratio = mgp.resb > 0. ? mgp.resa/mgp.resb : 0.;
  double mgpf_ratio = mgpf.resb > 0. ? mgpf.resa/mgpf.resb : 0.;
  fprintf (history_file, "%.12g,%.15g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%ld,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%d,%.12g,%.12g,%d,%.12g,%.12g,%.12g\n",
           t/t0, xcm, drift, umax, interface_T_mean, interface_T_rms,
           interface_gradT_tangent_rms, sigma_min, sigma_max, dt,
           grid->tn, mgp.resa, liquid_urms, liquid_umax, pcenter, pout,
           pcenter - pout, kapmean, kaprms, kapmin, kapmax, mgp.i,
           mgp.resb, mgp_ratio, mgpf.i, mgpf.resb, mgpf.resa, mgpf_ratio);
  fflush(history_file);
  if (i % 100 == 0)
    fprintf (stdout, "progress t/t0=%.6g dt=%.6g cells=%ld xcm=%.8g max_u=%.8g mass_drift=%.4g\n",
             t/t0, dt, grid->tn, xcm, umax, drift);
}

event adapt (i++)
{
  if (use_amr && (amr_mode != 2 || i == 0)) {
    if (amr_mode == 3)
      adapt_wavelet ({f,Tfield}, (double[]){1e-2,1e-5}, LEVEL);
    else
      adapt_wavelet ({f,u,Tfield}, (double[]){1e-2,1e-5,1e-5,1e-5}, LEVEL);
  }
}

event report (t = t_end)
{
  fprintf (stderr,
           "REPORT mode=%d PeT=%.8g dt_base=%.12g dt_factor=%.8g DT=%.12g tolerance=%.8g "
           "t_end_over_t0=%.8g mass_drift_max_pct=%.12g max_u=%.12g "
           "cells=%ld mgp_i=%d mgp_resa=%.12g\n",
           P1F_MODE, PeT, DT_BASE, dt_factor, DT, TOLERANCE, t/t0,
           100.*max_abs_mass_drift, last_max_velocity, grid->tn, mgp.i, mgp.resa);
  if (history_file)
    fclose(history_file), history_file = NULL;
}

event stop (t = t_end)
{
  return 1;
}
