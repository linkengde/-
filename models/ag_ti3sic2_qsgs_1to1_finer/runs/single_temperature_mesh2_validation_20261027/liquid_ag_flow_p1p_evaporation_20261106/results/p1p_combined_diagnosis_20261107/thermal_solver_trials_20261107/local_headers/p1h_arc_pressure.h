#ifndef P1H_ARC_PRESSURE_H
#define P1H_ARC_PRESSURE_H

#include <math.h>

typedef enum {
  P1H_PRESSURE_UNIFORM = 0,
  P1H_PRESSURE_PERIODIC_GAUSSIAN = 1
} P1HPressureProfile;

typedef struct {
  int enabled;
  double peak_Pa;
  double center_x_m, center_y_m;
  double radius_m;
  double start_s, ramp_s, end_s;
  double period_x_m, period_y_m;
  P1HPressureProfile profile;
} P1HArcPressure;

static inline double p1h_clamp01 (double x)
{
  return x < 0. ? 0. : x > 1. ? 1. : x;
}

static inline double p1h_smoothstep01 (double x)
{
  x = p1h_clamp01(x);
  return x*x*(3. - 2.*x);
}

static inline double p1h_heat_source_scale (int command_on, int ramping,
                                             double ramp_start_s,
                                             double ramp_s, double time_s)
{
  if (command_on) return 1.;
  if (!ramping || !(ramp_s > 0.)) return 0.;
  return 1. - p1h_smoothstep01((time_s-ramp_start_s)/ramp_s);
}

static inline double p1h_arc_pressure_time_factor (const P1HArcPressure *cfg,
                                                    double time_s)
{
  if (!cfg || !cfg->enabled || !(cfg->peak_Pa > 0.) ||
      time_s < cfg->start_s || time_s >= cfg->end_s)
    return 0.;
  double factor = cfg->ramp_s > 0. ?
    p1h_smoothstep01((time_s - cfg->start_s)/cfg->ramp_s) : 1.;
  if (cfg->ramp_s > 0. && isfinite(cfg->end_s))
    factor *= p1h_smoothstep01((cfg->end_s - time_s)/cfg->ramp_s);
  return factor;
}

static inline double p1h_periodic_delta (double x, double center, double period)
{
  double delta = x - center;
  if (period > 0.)
    delta -= period*nearbyint(delta/period);
  return delta;
}

static inline double p1h_arc_pressure_at (const P1HArcPressure *cfg,
                                           double x_m, double y_m,
                                           double time_s)
{
  double factor = p1h_arc_pressure_time_factor(cfg, time_s);
  if (factor == 0.)
    return 0.;
  double spatial = 1.;
  if (cfg->profile == P1H_PRESSURE_PERIODIC_GAUSSIAN) {
    if (!(cfg->radius_m > 0.))
      return 0.;
    double dx = p1h_periodic_delta(x_m, cfg->center_x_m, cfg->period_x_m);
    double dy = p1h_periodic_delta(y_m, cfg->center_y_m, cfg->period_y_m);
    spatial = exp(-(dx*dx + dy*dy)/(cfg->radius_m*cfg->radius_m));
  }
  return cfg->peak_Pa*factor*spatial;
}

/* f=1 is the sample and f=0 is gas, so the upper-surface grad(f) points
 * downward.  The positive pressure potential therefore creates a downward
 * force density p_eff*grad(f). */
static inline double p1h_arc_pressure_mask (double ag_fraction,
                                             double liquid_fraction,
                                             double normal_z)
{
  double top_facing = p1h_clamp01(-normal_z);
  return p1h_clamp01(ag_fraction)*p1h_clamp01(liquid_fraction)*top_facing;
}

static inline double p1h_arc_pressure_effective (const P1HArcPressure *cfg,
                                                  double x_m, double y_m,
                                                  double time_s,
                                                  double ag_fraction,
                                                  double liquid_fraction,
                                                  double normal_z)
{
  return p1h_arc_pressure_at(cfg, x_m, y_m, time_s)*
    p1h_arc_pressure_mask(ag_fraction, liquid_fraction, normal_z);
}

static inline double p1h_arc_pressure_force_density_z (double p_eff_Pa,
                                                        double grad_f_z_per_m)
{
  return p_eff_Pa*grad_f_z_per_m;
}

#endif
