/*
 * Static-circle check of the variable-sigma normal/tangential force split.
 * Uses Basilisk automatic height-function curvature; no exact-curvature path.
 * Dimensionless geometry only, with sigma(theta)=sigma0+sigma1*cos(theta).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "curvature.h"

scalar kappa[];
vector h[];
const double R = 0.5;
const double sigma0 = 1.;
const double sigma1 = 0.1;
int level_user = 4;

int main (int argc, char * argv[])
{
  if (argc != 2)
    return 2;
  level_user = atoi(argv[1]);
  if (level_user < 4 || level_user > 6)
    return 2;
  size (2.);
  origin (-1., -1.);
  init_grid (4*(1 << level_user));
  fraction (f, sq(R) - sq(x) - sq(y));
  boundary ({f});
  f.height = h;
  heights (f, h);
  cstats cs = curvature (f, kappa);

  double fnx = 0., fny = 0., ftx = 0., fty = 0.;
  double fn_abs = 0., ft_abs = 0., area = 0., ksum = 0., kn = 0.;
  foreach (reduction(+:fnx) reduction(+:fny) reduction(+:ftx)
           reduction(+:fty) reduction(+:fn_abs) reduction(+:ft_abs)
           reduction(+:area) reduction(+:ksum) reduction(+:kn)) {
    if (f[] > 0. && f[] < 1. && kappa[] < nodata) {
      coord n = interface_normal (point, f), pos;
      double alpha = plane_alpha (f[], n);
      double length = plane_area_center (n, alpha, &pos)*Delta;
      double xp = x + Delta*pos.x, yp = y + Delta*pos.y;
      double rr = sqrt(sq(xp) + sq(yp));
      double sig = sigma0 + sigma1*xp/rr;
      double gx = sigma1*sq(yp)/cube(rr);
      double gy = -sigma1*xp*yp/cube(rr);
      double nn = hypot(n.x,n.y), nx = n.x/nn, ny = n.y/nn;
      double ndotg = nx*gx + ny*gy;
      double tx = gx - ndotg*nx;
      double ty = gy - ndotg*ny;
      /* interface_normal points from liquid to gas; sigma*kappa*grad(f)
         therefore uses the inward normal -n for this positive-inside f. */
      double fn0 = -sig*kappa[]*nx*length;
      double fn1 = -sig*kappa[]*ny*length;
      double ft0 = tx*length;
      double ft1 = ty*length;
      fnx += fn0;
      fny += fn1;
      ftx += ft0;
      fty += ft1;
      fn_abs += hypot(fn0,fn1);
      ft_abs += hypot(ft0,ft1);
      area += length;
      ksum += kappa[]*length;
      kn += length;
    }
  }
  double close = hypot(fnx + ftx, fny + fty)/max(hypot(fnx,fny), hypot(ftx,fty));
  double normal_expected = -pi*sigma1;
  double tangential_expected = pi*sigma1;
  fprintf (stdout,
           "RESULT R_over_dx=%d N=%d sigma0=%.9g sigma1=%.9g normal_x=%.12g normal_y=%.12g normal_x_exact=%.12g normal_error_pct=%.12g tangential_x=%.12g tangential_y=%.12g tangential_x_exact=%.12g tangential_error_pct=%.12g force_closure_error_pct=%.12g normal_abs=%.12g tangential_abs=%.12g interface_length=%.12g mean_auto_curvature=%.12g curvature_calls_h=%d curvature_calls_fit=%d curvature_calls_avg=%d curvature_calls_centroid=%d\n",
           1 << level_user, 4*(1 << level_user), sigma0, sigma1,
           fnx, fny, normal_expected, 100.*fabs(fnx-normal_expected)/fabs(normal_expected),
           ftx, fty, tangential_expected, 100.*fabs(ftx-tangential_expected)/fabs(tangential_expected),
           100.*close, fn_abs, ft_abs, area, ksum/kn, cs.h, cs.f, cs.a, cs.c);
  return 0;
}
