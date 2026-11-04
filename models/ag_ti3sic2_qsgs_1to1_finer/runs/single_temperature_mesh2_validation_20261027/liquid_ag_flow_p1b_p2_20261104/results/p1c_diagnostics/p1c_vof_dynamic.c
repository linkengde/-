/*
 * P1b comparison: axisymmetric moving-drop Marangoni benchmark using the
 * project's automatic-curvature VOF path (two-phase.h + tension.h).
 * Dimensionless and independent of Ag material properties/geometry.
 * The imposed d(sigma)/dx is projected onto the local interface tangent.
 */
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

int LEVEL = 8;
int cells_per_radius = 16, force_mode = 0, use_amr = 1;
double timestep_factor = 1., tangent_scale = 1.;
const double R = 1. [1], NablaT = 1., Mu = 1., Rho = 1. [0];
const double Re = 0.066, Ca = 0.66;
const double Gamma_T = Re*sq(Mu)/(Rho*sq(R)*NablaT);
const double Gamma_0 = Gamma_T*R*NablaT/Ca;
const double t0 = Mu/(Gamma_T*NablaT);
const double Cdrop = 1., Cbulk = 1.;
double U_drop, u_drop, mass0, xcm0, tprev;
scalar deltaf[], nxf[], nyf[], kappa[];

int main (int argc, char * argv[])
{
  if (argc != 6) {
    fprintf (stderr, "usage: %s cells_per_radius force_mode(0=combined,1=normal,2=tangential) amr(0/1) dt_factor tangent_scale\n", argv[0]);
    return 2;
  }
  cells_per_radius = atoi(argv[1]);
  force_mode = atoi(argv[2]);
  use_amr = atoi(argv[3]);
  timestep_factor = atof(argv[4]);
  tangent_scale = atof(argv[5]);
  if (cells_per_radius < 8 || cells_per_radius > 32 ||
      force_mode < 0 || force_mode > 2 || (use_amr != 0 && use_amr != 1) || timestep_factor <= 0.)
    return 2;
  LEVEL = (int) (log(16.*cells_per_radius)/log(2.) + 0.5);
  N = 16*cells_per_radius;
  size (16.*R);
  origin (-L0/2.);
  rho1 = rho2 = Rho;
  mu1 = mu2 = Mu;
  f.sigma = Gamma_0 + fabs(Gamma_T)*R*NablaT;
  TOLERANCE = 1e-4 [*];
  /* Capillary stability bound for the requested uniform finest spacing. */
  DT = timestep_factor*sqrt(cube(R/cells_per_radius)/(pi*f.sigma));
  U_drop = -2./((2. + 3.*mu2/mu1)*(2. + Cdrop/Cbulk))*Gamma_T*R*NablaT/mu1;
  run();
}

static void interface_geometry()
{
  foreach() {
    double gx = (f[1] - f[-1])/(2.*Delta);
    double gy = (f[0,1] - f[0,-1])/(2.*Delta);
    double g = sqrt(sq(gx) + sq(gy));
    deltaf[] = g;
    nxf[] = g > 1e-30 ? gx/g : 0.;
    nyf[] = g > 1e-30 ? gy/g : 0.;
  }
  boundary ({deltaf,nxf,nyf});
}

event init (t = 0)
{
  /* f=0 is the spherical drop, matching the official test's (1-f) drop. */
  fraction (f, sq(x) + sq(y) - sq(R));
  boundary ({f});
  interface_geometry();
  foreach() {
    u.x[] = u.y[] = 0.;
    p[] = 0.;
  }
  boundary ((scalar *){u,p});
  mass0 = xcm0 = tprev = u_drop = 0.;
}

event acceleration (i++)
{
  interface_geometry();
  curvature (f, kappa, 1.);
  face vector av = a;
  foreach_face(x) {
    double ds = 0.5*(deltaf[] + deltaf[-1]);
    double tx = 0.5*((1. - sq(nxf[])) + (1. - sq(nxf[-1])));
    double kap = kappa[] < nodata && kappa[-1] < nodata ? 0.5*(kappa[] + kappa[-1]) : kappa[] < nodata ? kappa[] : kappa[-1];
    double grad_sigma = Gamma_T*NablaT;
    double dsigma = grad_sigma*x - fabs(grad_sigma)*R;
    if ((force_mode == 0 || force_mode == 1) && kap < nodata)
      av.x[] += alpha.x[]/(fm.x[] + SEPS)*dsigma*kap*(f[] - f[-1])/Delta;
    if (force_mode == 0 || force_mode == 2)
      av.x[] += alpha.x[]/(fm.x[] + SEPS)*tangent_scale*grad_sigma*tx*ds;
  }
  foreach_face(y) {
    double ds = 0.5*(deltaf[] + deltaf[0,-1]);
    double ty = -0.5*(nxf[]*nyf[] + nxf[0,-1]*nyf[0,-1]);
    double kap = kappa[] < nodata && kappa[0,-1] < nodata ? 0.5*(kappa[] + kappa[0,-1]) : kappa[] < nodata ? kappa[] : kappa[0,-1];
    double grad_sigma = Gamma_T*NablaT;
    double dsigma = grad_sigma*x - fabs(grad_sigma)*R;
    if ((force_mode == 0 || force_mode == 1) && kap < nodata)
      av.y[] += alpha.y[]/(fm.y[] + SEPS)*dsigma*kap*(f[] - f[0,-1])/Delta;
    if (force_mode == 0 || force_mode == 2)
      av.y[] += alpha.y[]/(fm.y[] + SEPS)*tangent_scale*grad_sigma*ty*ds;
  }
}

event logfile (i += 10)
{
  double xb = 0., vb = 0., sb = 0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb)) {
    double dvdrop = (1. - f[])*dv();
    xb += x*dvdrop;
    vb += u.x[]*dvdrop;
    sb += dvdrop;
  }
  if (t == 0.) {
    mass0 = sb;
    xcm0 = sb > 0. ? xb/sb : 0.;
    tprev = 0.;
    fprintf(stdout, "\nt mass_drift xb vb_over_Udrop tmid u_drop_over_Udrop dt cells mgp_i mgp_resa mode amr dtfac tangent_scale\n");
  }
  u_drop = t > tprev && sb > 0. ? (xb/sb - xcm0)/(t - tprev) : 0.;
  fprintf(stdout, "%g %.12g %.12g %.12g %g %.12g %.12g %ld %d %.12g %d %d %.8g %.8g\n",
          t/t0, mass0 > 0. ? (sb-mass0)/mass0 : 0.,
          sb > 0. ? xb/sb : 0., sb > 0. ? vb/sb/U_drop : 0.,
          (t+tprev)/2./t0, u_drop/U_drop, dt, grid->tn, mgp.i, mgp.resa, force_mode, use_amr, timestep_factor, tangent_scale);
  fflush(stdout);
  if (sb > 0.) xcm0 = xb/sb;
  tprev = t;
}

event report (t = 3.*t0)
{
  double U = -Gamma_T*R*NablaT/Mu;
  double umax = 0.;
  foreach(reduction(max:umax))
    umax = max(umax, norm(u));
  fprintf(stderr, "REPORT rdx=%d mode=%d amr=%d dt_factor=%.8g tangent_scale=%.8g measured_over_U=%.12g theory_over_U=%.12g rel_error_pct=%.12g max_velocity=%.12g mgp_i=%d mgp_resa=%.12g cells=%ld\n", cells_per_radius, force_mode, use_amr, timestep_factor, tangent_scale, u_drop/U, U_drop/U, 100.*fabs(u_drop/U_drop - 1.), umax, mgp.i, mgp.resa, grid->tn);
}

event stop (t = 3.*t0)
{
  return 1;
}

#if TREE
event adapt (i++)
{
  if (use_amr)
    adapt_wavelet ({f,u}, (double[]){1e-2,1e-5,1e-5}, LEVEL);
}
#endif
