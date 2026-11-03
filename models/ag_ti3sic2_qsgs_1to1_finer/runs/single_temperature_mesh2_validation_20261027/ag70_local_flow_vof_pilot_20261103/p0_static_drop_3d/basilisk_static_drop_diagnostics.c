/*
 * Dimensionless 3-D static-drop verification for P0.
 *
 * This is a numerical-method test only: rho = sigma = R = 1 and no
 * experimental Ag properties are implied.  `curvature_mode=exact` supplies
 * the analytical sphere curvature to isolate the pressure/force balance;
 * `curvature_mode=vof` uses Basilisk's VOF height-function curvature.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "vof.h"

attribute { double sigma; }
scalar f[], * interfaces = {f};
vector h[];

#include "iforce.h"

int max_level_user = 6;
int sessile = 0;
int exact_curvature = 0;
double theta_deg = 72.;
double dt_factor = 0.5;
double t_end = 0.5;
double mu0 = 1.;
double initial_volume = 0.;
double sphere_cz = 0.;
const double radius = 1.;
const double sigma0 = 1.;
const double box_size = 8.;

/* In Basilisk's 3-D boundary naming, z=Z0 is `back` (not `bottom`). */
h.t[back] = contact_angle(theta_deg*pi/180.);
h.r[back] = contact_angle(theta_deg*pi/180.);

int main (int argc, char * argv[])
{
  if (argc < 5) {
    fprintf (stderr, "usage: %s LEVEL sphere|sessile auto|exact DT_FACTOR [T_END] [MU]\n", argv[0]);
    return 2;
  }
  max_level_user = atoi(argv[1]);
  sessile = !strcmp(argv[2], "sessile");
  exact_curvature = !strcmp(argv[3], "exact");
  dt_factor = atof(argv[4]);
  if (argc > 5)
    t_end = atof(argv[5]);
  if (argc > 6)
    mu0 = atof(argv[6]);
  if (max_level_user < 6 || max_level_user > 9 || dt_factor <= 0. || t_end < 0. || mu0 < 0.)
    return 2;

  size (box_size);
  origin (-box_size/2., -box_size/2., sessile ? 0. : -box_size/2.);
  N = 4;
  /* 1e-4 is tighter than the <1% pressure-error acceptance target and
     avoids over-solving the viscous projection on the finest octree. */
  TOLERANCE = 1e-4;
  CFL = 0.5;
  stokes = true;
  DT = dt_factor*sqrt(cube(box_size/(1 << max_level_user))/(pi*sigma0));
  f.sigma = sigma0;
  f.height = h;
  const face vector muc[] = {mu0, mu0, mu0};
  mu = muc;
  run();
}

event init (i = 0)
{
  double ctheta = cos(theta_deg*pi/180.);
  double cz = sessile ? -radius*ctheta : sphere_cz;
  /* Refine geometrically around the implicit interface before fraction(). */
  refine (level < 3);
  for (int k = 0; k < max_level_user + 1; k++)
    refine (level < max_level_user &&
            fabs(sqrt(sq(x) + sq(y) + sq(z - cz)) - radius) < 3.*Delta);
  if (sessile) {
    fraction (f, sq(radius) - sq(x) - sq(y) - sq(z + radius*ctheta));
  }
  else
    fraction (f, sq(radius) - sq(x) - sq(y) - sq(z - sphere_cz));

  heights (f, h);
  initial_volume = statsf(f).sum;
  fprintf (stderr, "INIT kind=%s L0=%.12g origin=%.12g,%.12g,%.12g volume=%.12g domain_volume=%.12g cells=%ld DT=%.12g depth=%d\n",
           sessile ? "sessile72" : "sphere", L0, X0, Y0, Z0,
           initial_volume, statsf(f).volume, grid->tn, DT, grid->maxdepth);
  foreach() {
    u.x[] = u.y[] = u.z[] = 0.;
    /* Start from the Laplace-compatible pressure jump to avoid a startup impulse. */
    p[] = f[] > 0.5 ? 2.*sigma0/radius : 0.;
  }
  boundary ((scalar *){u, p});
}

/* Use the same interface-force discretisation as Basilisk's iforce.h. */
event surface_potential (i++, first)
{
  scalar phi = f.phi;
  if (!phi.i) {
    phi = new scalar;
    f.phi = phi;
  }
  if (exact_curvature) {
    foreach()
      phi[] = (f[] > 1e-6 && f[] < 1. - 1e-6) ? sigma0*2./radius : nodata;
    boundary ({phi});
  }
  else
    curvature (f, phi, sigma0, add = false);
}

event adapt (i++)
{
  scalar f1[];
  foreach()
    f1[] = f[];
  boundary ({f1});
  adapt_wavelet ({f1}, (double[]){1e-4}, minlevel = 3, maxlevel = max_level_user);
  heights (f, h);
}

static void diagnostics (void)
{
  double liquid_psum = 0., liquid_n = 0., gas_psum = 0., gas_n = 0.;
  double vmax = 0., liquid_u2 = 0., liquid_volume = 0.;
  double volume = statsf(f).sum;
  double initial = initial_volume;
  foreach(reduction(+:liquid_psum) reduction(+:liquid_n)
          reduction(+:gas_psum) reduction(+:gas_n)
          reduction(+:liquid_u2) reduction(+:liquid_volume)
          reduction(max:vmax)) {
    double ux = norm(u);
    if (ux > vmax) vmax = ux;
    double ff = clamp(f[], 0., 1.);
    liquid_u2 += ff*sq(ux)*dv();
    liquid_volume += ff*dv();
    double cz0 = sessile ? -radius*cos(theta_deg*pi/180.) : sphere_cz;
    double rr = sqrt(sq(x) + sq(y) + sq(z - cz0));
    if (rr < 0.72 && f[] > 0.99) {
      liquid_psum += p[];
      liquid_n += 1.;
    }
    if (rr > 1.25 && rr < 1.8 && f[] < 0.01) {
      gas_psum += p[];
      gas_n += 1.;
    }
  }
  double nfit = 0., sr = 0., sz = 0., srr = 0., szz = 0., srz = 0.;
  double sqsum = 0., qr = 0., qz = 0.;
  foreach(reduction(+:nfit) reduction(+:sr) reduction(+:sz)
          reduction(+:srr) reduction(+:szz) reduction(+:srz)
          reduction(+:sqsum) reduction(+:qr) reduction(+:qz)) {
    if (f[] > 0.01 && f[] < 0.99) {
      coord n = interface_normal (point, f), pos;
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &pos);
      double rxy = sqrt(sq(x + Delta*pos.x) + sq(y + Delta*pos.y));
      double zp = z + Delta*pos.z;
      double q = sq(rxy) + sq(zp);
      /* Algebraic circle fit: r^2 + z^2 = A*r + B*z + C. */
      nfit += 1.;
      sr += rxy;
      sz += zp;
      srr += sq(rxy);
      szz += sq(zp);
      srz += rxy*zp;
      sqsum += q;
      qr += q*rxy;
      qz += q*zp;
    }
  }
  double pjump = (liquid_n > 0. && gas_n > 0.) ? liquid_psum/liquid_n - gas_psum/gas_n : NAN;
  double volume_radius = cbrt(3.*volume/(4.*pi));
  double p_nominal = 2.*sigma0/radius;
  double p_volume = sessile ? NAN : 2.*sigma0/volume_radius;
  double pressure_error = 100.*(pjump - p_nominal)/p_nominal;
  double pressure_error_volume = sessile ? NAN : 100.*(pjump - p_volume)/p_volume;
  double analytical_volume = sessile ?
    pi*cube(radius)*(2. + cos(theta_deg*pi/180.))*sq(1. - cos(theta_deg*pi/180.))/3. :
    4.*pi*cube(radius)/3.;
  double initial_volume_error = 100.*(initial_volume - analytical_volume)/analytical_volume;
  double mass_drift = 100.*(volume - initial)/(initial > 0. ? initial : 1.);
  double theta_fit = NAN, fit_radius = NAN;
  if (sessile && nfit > 20.) {
    /* Centered least-squares fit of q=r^2+z^2 against r and z. */
    double covrr = srr - sr*sr/nfit;
    double covzz = szz - sz*sz/nfit;
    double covrz = srz - sr*sz/nfit;
    double covrq = qr - sr*sqsum/nfit;
    double covzq = qz - sz*sqsum/nfit;
    double det = covrr*covzz - covrz*covrz;
    if (fabs(det) > 1e-14) {
      double aa = (covrq*covzz - covzq*covrz)/det;
      double bb = (covzq*covrr - covrq*covrz)/det;
      double cc = (sqsum - aa*sr - bb*sz)/nfit;
      double zc = bb/2.;
      fit_radius = sqrt(fmax(0., sq(aa/2.) + sq(bb/2.) + cc));
      if (fit_radius > 0.)
        theta_fit = acos(fmax(-1., fmin(1., -zc/fit_radius)))*180./pi;
    }
  }
  double theta_err = sessile ? theta_fit - theta_deg : NAN;
  double urms_liquid = sqrt(liquid_u2/(liquid_volume > 0. ? liquid_volume : 1.));
  double ca = vmax*mu0/sigma0;
  double ca_rms = urms_liquid*mu0/sigma0;
  scalar kappa[];
  curvature (f, kappa);
  stats ks = statsf (kappa);
  double kmean = ks.volume > 0. ? ks.sum/ks.volume : NAN;
  double kmean_error = 100.*(kmean - 2.*sigma0/radius)/(2.*sigma0/radius);
  fprintf (stderr,
           "RESULT kind=%s curvature=%s level=%d R_over_dx=%d t=%.12g cells=%ld steps=%d dt=%.12g vmax=%.12g Ca=%.12g urms_liquid=%.12g Ca_rms=%.12g volume=%.12g initial_volume=%.12g initial_volume_error_pct=%.12g mass_drift_pct=%.12g p_jump=%.12g p_expected_nominal=%.12g pressure_error_pct=%.12g p_expected_volume=%.12g pressure_error_volume_pct=%.12g theta_deg=%.12g theta_error_deg=%.12g fitted_radius=%.12g kappa_mean=%.12g kappa_std=%.12g kappa_min=%.12g kappa_max=%.12g kappa_error_pct=%.12g\n",
           sessile ? "sessile72" : "sphere", exact_curvature ? "exact" : "vof", max_level_user,
           1 << (max_level_user - 3), t, grid->tn, iter, dt, vmax, ca, urms_liquid, ca_rms, volume, initial,
           initial_volume_error, mass_drift, pjump, p_nominal, pressure_error,
           p_volume, pressure_error_volume, theta_fit, theta_err, fit_radius,
           kmean, ks.stddev, ks.min, ks.max, kmean_error);
}

event timestep_diagnostics (i++, last)
{
  double vmax = 0., fmin = HUGE, fmax = -HUGE, pmin = HUGE, pmax = -HUGE;
  double liquid_u2 = 0., liquid_volume = 0.;
  int nonfinite_cells = 0;
  foreach (reduction(max:vmax) reduction(min:fmin) reduction(max:fmax)
           reduction(min:pmin) reduction(max:pmax)
           reduction(+:liquid_u2) reduction(+:liquid_volume)
           reduction(+:nonfinite_cells)) {
    double speed = norm(u), ff = clamp(f[], 0., 1.);
    if (speed > vmax) vmax = speed;
    if (f[] < fmin) fmin = f[];
    if (f[] > fmax) fmax = f[];
    if (p[] < pmin) pmin = p[];
    if (p[] > pmax) pmax = p[];
    if (!isfinite(f[]) || !isfinite(p[]) || !isfinite(u.x[]) ||
        !isfinite(u.y[]) || !isfinite(u.z[]))
      nonfinite_cells++;
    liquid_u2 += ff*sq(speed)*dv();
    liquid_volume += ff*dv();
  }
  double urms = sqrt(liquid_u2/(liquid_volume > 0. ? liquid_volume : 1.));
  double volume = statsf(f).sum;
  double mass_drift_pct = 100.*(volume - initial_volume)/(initial_volume > 0. ? initial_volume : 1.);
  double pressure_tolerance = TOLERANCE/sq(dt);
  double pressure_residual_ratio = mgp.resa/pressure_tolerance;
  fprintf (stderr, "STEP i=%d t=%.12g dt=%.12g pressure_iterations=%d pressure_resb=%.12g pressure_resa=%.12g pressure_tolerance=%.12g pressure_residual_ratio=%.12g pressure_rhs_sum=%.12g pressure_nrelax=%d vmax=%.12g urms_liquid=%.12g fmin=%.12g fmax=%.12g pmin=%.12g pmax=%.12g volume=%.12g mass_drift_pct=%.12g nonfinite_cells=%d cells=%ld\n",
           i, t, dt, mgp.i, mgp.resb, mgp.resa, pressure_tolerance, pressure_residual_ratio, mgp.sum, mgp.nrelax,
           vmax, urms, fmin, fmax, pmin, pmax, volume, mass_drift_pct, nonfinite_cells, grid->tn);
}

event report (t = t_end; t <= t_end)
{
  diagnostics();
  return 1;
}
