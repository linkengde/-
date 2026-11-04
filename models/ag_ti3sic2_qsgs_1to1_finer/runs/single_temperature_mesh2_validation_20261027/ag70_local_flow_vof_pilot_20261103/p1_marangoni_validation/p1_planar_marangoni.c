/*
 * Controlled 2-D planar-interface Marangoni benchmark.
 *
 * Dimensionless verification only; it does not use Ag properties or geometry.
 * A uniform prescribed d(sigma)/dx is applied as a CSF tangential force
 * d(sigma)/dx |grad(f)|.  For equal viscosities in two layers of height H,
 * no-slip walls, and a flat interface, the Stokes solution is
 * U_interface = (d(sigma)/dx) H/(2 mu).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

int level_user = 5;
double dt_factor = 0.5;
double t_end = 5.;
double sigma_ref = 0.01;
double tau_sigma = 0.;
double d_sigma_dT = 0.;
double grad_T = 0.01;
double H = 1.;
double Lx = 2.;
double initial_mass = 0.;
scalar delta_f[];

u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
uf.n[top] = 0.;
uf.n[bottom] = 0.;

static double row_mean (double yc)
{
  double sum = 0., n = 0.;
  foreach (reduction(+:sum) reduction(+:n))
    if (fabs(y - yc) < 0.1*Delta) {
      sum += u.x[];
      n += 1.;
    }
  return n > 0. ? sum/n : NAN;
}

static void compute_delta_f (void)
{
  foreach()
    delta_f[] = fabs((f[0,1] - f[0,-1])/(2.*Delta));
  boundary ({delta_f});
}

/* Put the prescribed interfacial shear into the centered solver's face force.
   The automatic curvature contribution continues to come from tension.h. */
event acceleration (i++)
{
  compute_delta_f();
  face vector av = a;
  foreach_face(x) {
    double ds = 0.5*(delta_f[] + delta_f[-1]);
    av.x[] += alpha.x[]/fm.x[]*tau_sigma*ds;
  }
}

event init (t = 0)
{
  fraction (f, -y);
  boundary ({f});
  compute_delta_f();
  initial_mass = statsf(f).sum;
  foreach() {
    u.x[] = u.y[] = 0.;
    p[] = 0.;
  }
  boundary ((scalar *){u,p});
  fprintf (stderr, "INIT level=%d H_over_dx=%d tau=%.12g sigma=%.12g mass=%.12g cells=%ld\n",
           level_user, 1 << level_user, tau_sigma, sigma_ref, initial_mass,
           grid->tn);
}

event sample (t += 1.)
{
  compute_delta_f();
  double un = 0., wr = 0., u2 = 0., umax = 0., mass = 0., force_weight = 0.;
  foreach (reduction(+:un) reduction(+:wr) reduction(+:u2)
           reduction(max:umax) reduction(+:mass) reduction(+:force_weight)) {
    double w = delta_f[]*dv();
    un += u.x[]*w;
    u2 += sq(u.x[])*w;
    wr += w;
    force_weight += w;
    umax = max(umax, norm(u));
    mass += f[]*dv();
  }
  double umean = wr > 0. ? un/wr : NAN;
  double urms = wr > 0. ? sqrt(u2/wr) : NAN;
  fprintf (stderr,
           "SAMPLE level=%d t=%.9g tau=%.9g u_interface=%.12g u_rms_interface=%.12g umax=%.12g mass_drift_pct=%.12g force_weight_per_length=%.12g mgp_i=%d mgp_resa=%.9g\n",
           level_user, t, tau_sigma, umean, urms, umax,
           100.*(mass-initial_mass)/initial_mass, force_weight/Lx,
           mgp.i, mgp.resa);
}

event final_report (t = t_end; t <= t_end)
{
  compute_delta_f();
  double mass = statsf(f).sum;
  double umax = 0., u_num = 0., u2 = 0., wsum = 0., force_weight = 0.;
  foreach (reduction(max:umax) reduction(+:u_num) reduction(+:u2) reduction(+:wsum)
           reduction(+:force_weight)) {
    double w = delta_f[]*dv();
    u_num += u.x[]*w;
    u2 += sq(u.x[])*w;
    wsum += w;
    force_weight += w;
    umax = max(umax, norm(u));
  }
  double u_interface = wsum > 0. ? u_num/wsum : NAN;
  double Delta0 = H/(1 << level_user);
  double y0 = 0.5*Delta0;
  double up0 = row_mean(y0), up1 = row_mean(3.*y0);
  double um0 = row_mean(-y0), um1 = row_mean(-3.*y0);
  double shear_top = (up1-up0)/Delta0;
  double shear_bottom = (um0-um1)/Delta0;
  double viscous_jump = mu1*(shear_top-shear_bottom);
  double force_per_length = tau_sigma*force_weight/Lx;
  double balance_error = fabs(tau_sigma + viscous_jump)/max(fabs(tau_sigma),1e-30);
  double u_exact = tau_sigma*H/(2.*mu1);
  double rel_error = tau_sigma != 0. ? 100.*fabs(u_interface-u_exact)/fabs(u_exact) : NAN;
  double mass_drift = 100.*(mass-initial_mass)/initial_mass;
  fprintf (stderr,
           "RESULT level=%d H_over_dx=%d N=%d mode=%s dt_factor=%.9g t=%.9g tau=%.12g d_sigma_dT=%.12g grad_T=%.12g u_interface=%.12g u_exact=%.12g velocity_error_pct=%.12g u_rms_interface=%.12g umax=%.12g mass_drift_pct=%.12g force_per_length=%.12g force_expected=%.12g viscous_stress_jump=%.12g stress_balance_error_pct=%.12g shear_top=%.12g shear_bottom=%.12g mgp_i=%d mgp_resa=%.12g cells=%ld\n",
           level_user, 1 << level_user, 2*(1 << level_user),
           tau_sigma == 0. ? "zero_gradient" : tau_sigma < 0. ? "negative_dsigma_dT" : "positive_dsigma_dT",
           dt_factor, t, tau_sigma, d_sigma_dT, grad_T, u_interface, u_exact,
           rel_error, wsum > 0. ? sqrt(u2/wsum) : NAN, umax, mass_drift,
           force_per_length, tau_sigma, viscous_jump, 100.*balance_error,
           shear_top, shear_bottom, mgp.i, mgp.resa, grid->tn);
  return 1;
}

int main (int argc, char * argv[])
{
  if (argc < 4) {
    fprintf (stderr, "usage: %s LEVEL zero|neg|pos|agred|agoxygen DT_FACTOR [T_END]\n", argv[0]);
    return 2;
  }
  level_user = atoi(argv[1]);
  dt_factor = atof(argv[3]);
  if (argc > 4)
    t_end = atof(argv[4]);
  if (level_user < 4 || level_user > 7 || dt_factor <= 0. || t_end <= 0.)
    return 2;
  const char * mode = argv[2];
  if (!strcmp(mode, "neg"))
    d_sigma_dT = -0.5;
  else if (!strcmp(mode, "pos"))
    d_sigma_dT = 0.5;
  else if (!strcmp(mode, "agred"))
    /* Morohoshi et al. (2008): -0.298 mN/(m K), reducing Ar/He + 5% H2. */
    d_sigma_dT = 0.5*(-0.298/0.546);
  else if (!strcmp(mode, "agoxygen"))
    /* Morohoshi et al. (2008): +0.546 mN/(m K), pO2 = 100 Pa. */
    d_sigma_dT = 0.5;
  else if (strcmp(mode, "zero"))
    return 2;
  tau_sigma = d_sigma_dT*grad_T;

  size (Lx);
  origin (0., -H);
  init_grid (2*(1 << level_user));
  periodic (right);
  rho1 = rho2 = 1.;
  mu1 = mu2 = 1.;
  f.sigma = sigma_ref;
  TOLERANCE = 1e-7;
  CFL = 0.5;
  stokes = true;
  DT = dt_factor*sqrt(cube(1./(1 << level_user))/(pi*sigma_ref));
  run();
  return 0;
}
