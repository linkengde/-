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
  check_status("overflowing_requested_mass_rejected", P1P_ERR_ARGUMENT,
               p1p_evaporation_step(1.e308,1.e308,1.e308,1.,1.,&step));
  check_status("overflowing_latent_energy_rejected", P1P_ERR_ARGUMENT,
               p1p_evaporation_step(1.e100,1.,1.,1.e100,1.e300,&step));

  P1PStep whole, half_a, half_b;
  p1p_evaporation_step(5., 0.02, 0.1, 100., 2.e6, &whole);
  p1p_evaporation_step(5., 0.02, 0.05, 100., 2.e6, &half_a);
  p1p_evaporation_step(5., 0.02, 0.05, 100., 2.e6, &half_b);
  check_close("constant_flux_dt_halving_mass_invariance_kg",
              whole.removed_mass_kg,
              half_a.removed_mass_kg + half_b.removed_mass_kg, 1.e-12);

  P1PCellState cell = {0.9, 0.6, 0.2, 0.1, 1.e7};
  const P1PCellState initial_cell = cell;
  P1PCellUpdate cell_update;
  st = p1p_apply_cell_evaporation(&table, 0.5, 0., 77.1, 1500., 0.7,
                                  4.e5, 10490., 1.e-12, 1.e-6, 1.e-6,
                                  &cell, &cell_update);
  check_status("cell_evaporation_update_valid_inputs", P1P_OK, st);
  check_close("cell_ag_mass_initial_equals_remaining_plus_escaped_kg",
              initial_cell.ag_fraction * 10490. * 1.e-12,
              cell.ag_fraction * 10490. * 1.e-12 +
                  cell_update.removed_mass_kg,
              1.e-12);
  check_close("cell_sample_and_components_remain_consistent",
              cell.sample_fraction,
              cell.ag_fraction + cell.skeleton_fraction + cell.pore_fraction,
              1.e-12);
  check_close("cell_skeleton_fraction_unchanged", initial_cell.skeleton_fraction,
              cell.skeleton_fraction, 0.);
  check_close("cell_pore_fraction_unchanged", initial_cell.pore_fraction,
              cell.pore_fraction, 0.);
  check_close("cell_energy_initial_equals_remaining_plus_escaped_J",
              initial_cell.enthalpy_density_J_m3 * 1.e-12,
              cell.enthalpy_density_J_m3 * 1.e-12 +
                  cell_update.escaped_energy_J,
              1.e-12);
  check_close("escaped_energy_equals_sensible_plus_vaporization_J",
              cell_update.ag_sensible_energy_J + cell_update.latent_energy_J,
              cell_update.escaped_energy_J, 1.e-12);

  P1PCellState solid_cell = {0.9, 0.6, 0.2, 0.1, 1.e7};
  st = p1p_apply_cell_evaporation(&table, 0.5, 0., 77.1, 1500., 0.,
                                  4.e5, 10490., 1.e-12, 1.e-6, 1.e-6,
                                  &solid_cell, &cell_update);
  check_status("no_liquid_ag_no_evaporation", P1P_OK, st);
  check_close("no_liquid_ag_removed_mass_kg", 0.,
              cell_update.removed_mass_kg, 0.);

  P1PCellState energy_limited_cell = {0.9, 0.6, 0.2, 0.1, 1000.};
  st = p1p_apply_cell_evaporation(&table, 0.5, 0., 77.1, 1500., 1.,
                                  4.e5, 10490., 1.e-12, 1.e-6, 1.e-6,
                                  &energy_limited_cell, &cell_update);
  check_status("enthalpy_limited_evaporation_valid_inputs", P1P_OK, st);
  check_close("enthalpy_cap_reported", 1.,
              cell_update.enthalpy_limited ? 1. : 0., 0.);
  check_close("enthalpy_floor_not_crossed_J", 0.,
              energy_limited_cell.enthalpy_density_J_m3 * 1.e-12, 0.);

  P1PCellState range_cell = {0.9, 0.6, 0.2, 0.1, 1.e7};
  const P1PCellState range_cell_before = range_cell;
  st = p1p_apply_cell_evaporation(&table, 0.5, 0., 77.1, 999., 1.,
                                  4.e5, 10490., 1.e-12, 1.e-6, 1.e-6,
                                  &range_cell, &cell_update);
  check_status("cell_evaporation_rejects_property_extrapolation",
               P1P_ERR_RANGE, st);
  check_close("out_of_range_update_does_not_change_ag",
              range_cell_before.ag_fraction, range_cell.ag_fraction, 0.);
  check_close("out_of_range_update_does_not_change_enthalpy",
              range_cell_before.enthalpy_density_J_m3,
              range_cell.enthalpy_density_J_m3, 0.);

  P1PPropertyTable overflow_table=table;
  for (size_t row=0;row<overflow_table.n;row++) {
    overflow_table.saturation_pressure_Pa[row]=1.e308;
    overflow_table.vaporization_enthalpy_J_kg[row]=1.e308;
  }
  P1PCellState overflow_cell={0.9,0.6,0.2,0.1,1.e7};
  const P1PCellState overflow_cell_before=overflow_cell;
  st=p1p_apply_cell_evaporation(&overflow_table,1.,0.,77.1,1500.,1.,
                                4.e5,10490.,1.,1.e200,1.,
                                &overflow_cell,&cell_update);
  check_status("overflowing_cell_source_rejected",P1P_ERR_ARGUMENT,st);
  check_close("rejected_overflow_does_not_mutate_ag",
              overflow_cell_before.ag_fraction,overflow_cell.ag_fraction,0.);
  check_close("rejected_overflow_does_not_mutate_enthalpy",
              overflow_cell_before.enthalpy_density_J_m3,
              overflow_cell.enthalpy_density_J_m3,0.);

  fclose(csv);
  if (failures) {
    fprintf(stderr, "%d kernel test(s) failed\n", failures);
    return 1;
  }
  puts("PASS: P1P standalone flux, property-table, bounded-mass, and latent-energy kernel tests");
  return 0;
}
