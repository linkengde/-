/*
 * Reconstruct the actual P1b axisymmetric face-CSF normal and tangential
 * force split on a spherical interface and compare their volume-integrated
 * resultants with the analytic closed-surface balance.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

int RDX = 8;
const double R = 1., G = 0.066, sigma0 = 0.1;
scalar deltaf[], nxf[], nyf[], kappa[];
face vector fn[], ft[];

static void interface_geometry (void)
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

int main (int argc, char * argv[])
{
  if (argc != 2) return 2;
  RDX = atoi(argv[1]);
  if (RDX != 8 && RDX != 16 && RDX != 32) return 2;
  size (16.);
  origin (-8., 0.);
  init_grid (16*RDX);
  if (is_constant(cm)) {
    scalar * l = list_copy (all);
    cm = new scalar;
    free (all);
    all = list_concat ({cm}, l);
    free (l);
  }
  scalar cmv = cm;
  foreach() cmv[] = y;
  if (is_constant(fm.x)) {
    scalar * l = list_copy (all);
    fm = new face vector;
    free (all);
    all = list_concat ((scalar *){fm}, l);
    free (l);
  }
  face vector fmv = fm;
  foreach_face() fmv.x[] = max(y, 1./HUGE);
  fraction (f, sq(x) + sq(y) - sq(R));
  boundary ({f});
  interface_geometry();
  curvature (f, kappa, 1.);
  foreach_face() fn.x[] = ft.x[] = 0.;
  foreach_face(x) {
    double ds = 0.5*(deltaf[] + deltaf[-1]);
    double tx = 0.5*((1. - sq(nxf[])) + (1. - sq(nxf[-1])));
    double kap = kappa[] < nodata && kappa[-1] < nodata ? 0.5*(kappa[] + kappa[-1]) : kappa[] < nodata ? kappa[] : kappa[-1];
    double dsigma = G*x - fabs(G)*R;
    if (kap < nodata)
      fn.x[] += alpha.x[]/(fm.x[] + SEPS)*dsigma*kap*(f[] - f[-1])/Delta;
    ft.x[] += alpha.x[]/(fm.x[] + SEPS)*G*tx*ds;
  }
  foreach_face(y) {
    double ds = 0.5*(deltaf[] + deltaf[0,-1]);
    double ty = -0.5*(nxf[]*nyf[] + nxf[0,-1]*nyf[0,-1]);
    double kap = kappa[] < nodata && kappa[0,-1] < nodata ? 0.5*(kappa[] + kappa[0,-1]) : kappa[] < nodata ? kappa[] : kappa[0,-1];
    double dsigma = G*x - fabs(G)*R;
    if (kap < nodata)
      fn.y[] += alpha.y[]/(fm.y[] + SEPS)*dsigma*kap*(f[] - f[0,-1])/Delta;
    ft.y[] += alpha.y[]/(fm.y[] + SEPS)*G*ty*ds;
  }
  double fnx=0., fny=0., ftx=0., fty=0., maxk=0., nk=0.;
  foreach (reduction(+:fnx) reduction(+:fny) reduction(+:ftx) reduction(+:fty) reduction(max:maxk) reduction(+:nk)) {
    fnx += 0.5*(fn.x[]*fm.x[] + fn.x[1]*fm.x[1])*dv();
    fny += 0.5*(fn.y[]*fm.y[] + fn.y[0,1]*fm.y[0,1])*dv();
    ftx += 0.5*(ft.x[]*fm.x[] + ft.x[1]*fm.x[1])*dv();
    fty += 0.5*(ft.y[]*fm.y[] + ft.y[0,1]*fm.y[0,1])*dv();
    if (kappa[] < nodata) { maxk = max(maxk, fabs(kappa[])); nk++; }
  }
  /* 2*pi is omitted in Basilisk's axisymmetric metric; exact F/(2*pi)=4/3 G R^2. */
  double exact = 4.*G*sq(R)/3.;
  double closure = fabs(fnx + ftx)/max(fabs(fnx), fabs(ftx));
  fprintf (stdout, "R_over_dx,N,normal_x_over_2pi,normal_r_over_2pi,tangent_x_over_2pi,tangent_r_over_2pi,analytic_abs_axial_each_over_2pi,normal_axial_error_pct,tangent_axial_error_pct,axial_force_closure_pct,max_abs_kappa,curvature_cells\n");
  fprintf (stdout, "%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.9g,%.9g,%.9g,%.12g,%.0f\n", RDX, 16*RDX, fnx, fny, ftx, fty, exact, 100.*fabs(fabs(fnx)-exact)/exact, 100.*fabs(fabs(ftx)-exact)/exact, 100.*closure, maxk, nk);
  return 0;
}
