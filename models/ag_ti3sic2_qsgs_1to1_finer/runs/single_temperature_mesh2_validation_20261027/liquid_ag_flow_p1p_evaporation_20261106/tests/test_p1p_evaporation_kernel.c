#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../p1p_evaporation.h"

static int failures = 0;
static FILE *csv;

static void check_close(const char *name, double expected, double actual,
                        double rel_tol)
{
  double err = fabs(actual - expected) / fmax(fabs(expected), 1.e-300);
  int ok = isfinite(actual) && err <= rel_tol;
  fprintf(csv, "%s,%.17g,%.17g,%.17g,%s\n", name, expected, actual, err,
          ok ? "PASS" : "FAIL");
  if (!ok) {
    fprintf(stderr, "FAIL %s expected=%.17g actual=%.17g relerr=%.6g\n",
            name, expected, actual, err);
    failures++;
  }
}

static void check_status(const char *name, P1PStatus expected,
                         P1PStatus actual)
{
  int ok = actual == expected;
  fprintf(csv, "%s,%d,%d,0,%s\n", name, (int)expected, (int)actual,
          ok ? "PASS" : "FAIL");
  if (!ok) {
    fprintf(stderr, "FAIL %s expected_status=%d actual_status=%d\n",
            name, (int)expected, (int)actual);
    failures++;
  }
}

int main(int argc, char **argv)
{
  if (argc != 3) {
    fprintf(stderr, "usage: %s synthetic_property_table.csv output.csv\n", argv[0]);
    return 2;
  }
  csv = fopen(argv[2], "w");
  if (!csv) {
    perror(argv[2]);
    return 2;
  }
  fputs("test,expected,actual,relative_error,status\n", csv);

  P1PPropertyTable table;
  check_status("synthetic_table_load", P1P_OK,
               p1p_load_property_table(argv[1], &table));
  double psat = NAN, lv = NAN;
  check_status("synthetic_table_log_pressure_lookup", P1P_OK,
               p1p_property_lookup(&table, 1500., &psat, &lv));
  check_close("synthetic_geometric_midpoint_pressure_Pa", sqrt(1000.*10000.),
              psat, 1.e-12);
  check_close("synthetic_linear_midpoint_Lv_J_kg", 1900000., lv, 1.e-12);
  check_status("no_property_extrapolation", P1P_ERR_RANGE,
               p1p_property_lookup(&table, 999., &psat, &lv));

  const double alpha = 0.5, p_sat = 10000., p_v = 1000.;
  const double t = 1500., r_specific = 77.1;
  double flux = NAN;
  P1PStatus st = p1p_hk_mass_flux(alpha, p_sat, p_v, t, r_specific, &flux);
  check_status("hertz_knudsen_valid_inputs", P1P_OK, st);
  check_close("hertz_knudsen_analytic_flux_kg_m2_s",
              alpha * (p_sat-p_v) / sqrt(2.*P1P_PI*r_specific*t), flux, 1.e-12);
  p1p_hk_mass_flux(1., 5000., 5000., t, r_specific, &flux);
  check_close("zero_pressure_difference_flux", 0., flux, 0.);
  p1p_hk_mass_flux(1., 5000., 6000., t, r_specific, &flux);
  check_close("negative_pressure_difference_clipped_to_zero", 0., flux, 0.);
  check_status("invalid_evaporation_coefficient_rejected", P1P_ERR_ARGUMENT,
               p1p_hk_mass_flux(1.01, p_sat, p_v, t, r_specific, &flux));

  P1PStep step;
  st = p1p_evaporation_step(2., 0.5, 2., 0.8, 2.e6, &step);
  check_status("mass_limited_step_valid_inputs", P1P_OK, st);
  check_close("requested_mass_kg", 2., step.requested_mass_kg, 1.e-12);
  check_close("removed_mass_capped_to_available_kg", 0.8, step.removed_mass_kg, 1.e-12);
  check_close("latent_energy_J", 1.6e6, step.latent_energy_J, 1.e-12);
  check_close("mass_balance_initial_equals_remaining_plus_escaped_kg",
              0.8, (0.8-step.removed_mass_kg)+step.removed_mass_kg, 1.e-12);
  check_close("mass_limit_flag", 1., step.mass_limited ? 1. : 0., 0.);

  P1PStep whole, half_a, half_b;
  p1p_evaporation_step(5., 0.02, 0.1, 100., 2.e6, &whole);
  p1p_evaporation_step(5., 0.02, 0.05, 100., 2.e6, &half_a);
  p1p_evaporation_step(5., 0.02, 0.05, 100., 2.e6, &half_b);
  check_close("constant_flux_dt_halving_mass_invariance_kg",
              whole.removed_mass_kg,
              half_a.removed_mass_kg + half_b.removed_mass_kg, 1.e-12);

  fclose(csv);
  if (failures) {
    fprintf(stderr, "%d kernel test(s) failed\n", failures);
    return 1;
  }
  puts("PASS: P1P standalone flux, property-table, bounded-mass, and latent-energy kernel tests");
  return 0;
}
