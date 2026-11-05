#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "p1h_arc_pressure.h"

int main (void)
{
  P1HArcPressure cfg = {
    .enabled = 1, .peak_Pa = 100.,
    .center_x_m = 0., .center_y_m = 0., .radius_m = 1.,
    .start_s = 0., .ramp_s = 1., .end_s = INFINITY,
    .period_x_m = 10., .period_y_m = 10.,
    .profile = P1H_PRESSURE_UNIFORM
  };
  assert(p1h_arc_pressure_at(&cfg, 2., 3., 0.) == 0.);
  assert(fabs(p1h_arc_pressure_at(&cfg, 2., 3., 0.5) - 50.) < 1.e-12);
  assert(fabs(p1h_arc_pressure_at(&cfg, 2., 3., 1.) - 100.) < 1.e-12);
  assert(p1h_arc_pressure_mask(1., 1., -1.) == 1.);
  assert(p1h_arc_pressure_mask(0.5, 1., -1.) == 0.5);
  assert(p1h_arc_pressure_mask(1., 0., -1.) == 0.);
  assert(p1h_arc_pressure_mask(1., 1., 1.) == 0.);
  assert(p1h_arc_pressure_mask(1., 1., 0.) == 0.);
  assert(p1h_arc_pressure_force_density_z(100., -2.) == -200.);

  /* A source cut and a pressure cut using the same off-ramp stay synchronized. */
  cfg.profile = P1H_PRESSURE_UNIFORM;
  cfg.start_s = 0.; cfg.ramp_s = 2.; cfg.end_s = 4.; cfg.peak_Pa = 100.;
  for (int n = 0; n <= 8; n++) {
    double time = 2. + 0.25*n;
    double source_scale = p1h_heat_source_scale(0,1,2.,2.,time);
    double pressure_scale = p1h_arc_pressure_time_factor(&cfg,time);
    assert(fabs(source_scale-pressure_scale) < 1.e-12);
  }

  cfg.profile = P1H_PRESSURE_PERIODIC_GAUSSIAN;
  cfg.ramp_s = 0.;
  cfg.radius_m = 1.;
  double p0 = p1h_arc_pressure_at(&cfg, 0.2, 0.3, 0.5);
  double p1 = p1h_arc_pressure_at(&cfg, 10.2, 0.3, 0.5);
  assert(fabs(p0 - p1) < 1.e-12); /* periodic seam continuity */
  assert(p1h_arc_pressure_at(&cfg, 5.2, 0.3, 0.5) < p0);

  /* A planar upper VOF interface has grad(f).z < 0 because f=1 is sample. */
  const int n = 32;
  const double dx = 1.e-9, p = 100.;
  double integrated_force_z = 0.;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      double normal_z = -1.;
      double mask = p1h_arc_pressure_mask(1., 1., normal_z);
      integrated_force_z += p1h_arc_pressure_force_density_z(p*mask, -1./dx)*dx*dx*dx;
    }
  assert(integrated_force_z < 0.);
  printf("P1H pressure module tests passed; planar upper-surface Fz=%.9g N\n",
         integrated_force_z);
  return 0;
}
