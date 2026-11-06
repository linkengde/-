#ifndef P1P_EVAPORATION_H
#define P1P_EVAPORATION_H

/*
 * Isolated P1P Hertz-Knudsen evaporation kernel.
 *
 * This header deliberately contains no recoil-pressure term. It accepts a
 * sourced saturation-pressure/latent-heat table and reports a bounded local
 * mass and latent-energy removal for a caller-supplied interface area.
 */
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define P1P_MAX_PROPERTY_ROWS 256
#define P1P_PI 3.14159265358979323846264338327950288

typedef enum {
  P1P_OK = 0,
  P1P_ERR_ARGUMENT = 1,
  P1P_ERR_TABLE = 2,
  P1P_ERR_RANGE = 3
} P1PStatus;

typedef struct {
  size_t n;
  double temperature_K[P1P_MAX_PROPERTY_ROWS];
  double saturation_pressure_Pa[P1P_MAX_PROPERTY_ROWS];
  double vaporization_enthalpy_J_kg[P1P_MAX_PROPERTY_ROWS];
  /* NaN means alpha_ev is unavailable at this row. Finite alpha anchors are
   * interpolated only inside their source-supported temperature interval. */
  double alpha_ev[P1P_MAX_PROPERTY_ROWS];
} P1PPropertyTable;

typedef struct {
  double mass_flux_kg_m2_s;
  double requested_mass_kg;
  double removed_mass_kg;
  double latent_energy_J;
  int mass_limited;
} P1PStep;

/* Cell-averaged conservative state used by the solver source hook.  The
 * fractions are volume fractions in the complete cell, with sample_fraction
 * = Ag + skeleton + unresolved pore.  Enthalpy is a bulk-cell density in
 * J/m^3, referenced to 300 K by the P1M material model. */
typedef struct {
  double sample_fraction;
  double ag_fraction;
  double skeleton_fraction;
  double pore_fraction;
  double enthalpy_density_J_m3;
} P1PCellState;

typedef struct {
  double mass_flux_kg_m2_s;
  double requested_mass_kg;
  double removed_mass_kg;
  double ag_sensible_energy_J;
  double latent_energy_J;
  double escaped_energy_J;
  double removed_volume_fraction;
  double alpha_ev_used;
  int alpha_from_property_table;
  int liquid_mass_limited;
  int enthalpy_limited;
} P1PCellUpdate;

static inline int p1p_parse_csv_number(char *field, double *value,
                                       int allow_blank)
{
  while (*field == ' ' || *field == '\t' || *field == '\r' || *field == '\n')
    field++;
  char *end = field + strlen(field);
  while (end > field && (end[-1] == ' ' || end[-1] == '\t' ||
                         end[-1] == '\r' || end[-1] == '\n'))
    *--end = '\0';
  if (!*field) {
    if (!allow_blank) return 0;
    *value = NAN;
    return 1;
  }
  char *tail = NULL;
  double parsed = strtod(field, &tail);
  if (tail == field || *tail != '\0' || !isfinite(parsed)) return 0;
  *value = parsed;
  return 1;
}

/* CSV contracts: legacy temperature_K,p_sat_Pa,L_v_J_kg, or the same header
 * with a fourth alpha_ev column. Empty alpha cells mean no value is available. */
static inline P1PStatus p1p_load_property_table(const char *path,
                                                P1PPropertyTable *table)
{
  if (!path || !table)
    return P1P_ERR_ARGUMENT;
  FILE *fp = fopen(path, "r");
  if (!fp)
    return P1P_ERR_TABLE;
  memset(table, 0, sizeof(*table));
  char line[512];
  if (!fgets(line, sizeof(line), fp)) {
    fclose(fp);
    return P1P_ERR_TABLE;
  }
  if ((unsigned char)line[0] == 0xef && (unsigned char)line[1] == 0xbb &&
      (unsigned char)line[2] == 0xbf)
    memmove(line, line + 3, strlen(line + 3) + 1);
  char *newline = strpbrk(line, "\r\n");
  if (newline) *newline = '\0';
  static const char *header3 = "temperature_K,p_sat_Pa,L_v_J_kg";
  static const char *header4 = "temperature_K,p_sat_Pa,L_v_J_kg,alpha_ev";
  int has_alpha_column = strcmp(line, header4) == 0;
  if (!has_alpha_column && strcmp(line, header3) != 0) {
    fclose(fp);
    return P1P_ERR_TABLE;
  }
  while (fgets(line, sizeof(line), fp)) {
    if (line[0] == '\n' || line[0] == '\r' || line[0] == '#')
      continue;
    if (table->n >= P1P_MAX_PROPERTY_ROWS) {
      fclose(fp);
      return P1P_ERR_TABLE;
    }
    char *fields[4] = {NULL, NULL, NULL, NULL};
    fields[0] = line;
    int field_count = 1;
    for (char *cursor = line; *cursor && field_count < 4; cursor++) {
      if (*cursor == ',') {
        *cursor = '\0';
        fields[field_count++] = cursor + 1;
      }
    }
    double t, p, lv, alpha = NAN;
    if (field_count != (has_alpha_column ? 4 : 3) ||
        !p1p_parse_csv_number(fields[0], &t, 0) ||
        !p1p_parse_csv_number(fields[1], &p, 0) ||
        !p1p_parse_csv_number(fields[2], &lv, 0) ||
        (has_alpha_column && !p1p_parse_csv_number(fields[3], &alpha, 1)) ||
        !isfinite(t) || !isfinite(p) || !isfinite(lv) ||
        (has_alpha_column && isfinite(alpha) && (alpha < 0. || alpha > 1.)) ||
        t <= 0. || p <= 0. || lv <= 0. ||
        (table->n && t <= table->temperature_K[table->n - 1])) {
      fclose(fp);
      return P1P_ERR_TABLE;
    }
    size_t i = table->n++;
    table->temperature_K[i] = t;
    table->saturation_pressure_Pa[i] = p;
    table->vaporization_enthalpy_J_kg[i] = lv;
    table->alpha_ev[i] = alpha;
  }
  fclose(fp);
  return table->n >= 2 ? P1P_OK : P1P_ERR_TABLE;
}

/* Interpolate alpha_ev only between finite source anchors; never extrapolate.
 * A separately declared fallback may be used for sensitivity cases. */
static inline P1PStatus p1p_alpha_lookup(const P1PPropertyTable *table,
                                         double temperature_K,
                                         double *alpha_ev)
{
  if (!table || !alpha_ev || !isfinite(temperature_K) ||
      temperature_K <= 0. || table->n < 1)
    return P1P_ERR_ARGUMENT;
  size_t first = table->n, previous = table->n;
  for (size_t i = 0; i < table->n; i++) {
    if (!isfinite(table->alpha_ev[i])) continue;
    if (first == table->n) first = i;
    if (temperature_K == table->temperature_K[i]) {
      *alpha_ev = table->alpha_ev[i];
      return P1P_OK;
    }
    if (table->temperature_K[i] < temperature_K) previous = i;
    if (table->temperature_K[i] > temperature_K) {
      if (previous == table->n) return P1P_ERR_RANGE;
      const double t0 = table->temperature_K[previous];
      const double t1 = table->temperature_K[i];
      const double f = (temperature_K - t0)/(t1 - t0);
      *alpha_ev = table->alpha_ev[previous] +
          f*(table->alpha_ev[i] - table->alpha_ev[previous]);
      return P1P_OK;
    }
  }
  if (first == table->n || temperature_K < table->temperature_K[first])
    return P1P_ERR_RANGE;
  return P1P_ERR_RANGE;
}

/*
 * Interpolate ln(p_sat) linearly in T and L_v linearly in T. No extrapolation
 * is permitted; outside-table temperatures fail closed.
 */
static inline P1PStatus p1p_property_lookup(const P1PPropertyTable *table,
                                            double temperature_K,
                                            double *p_sat_Pa,
                                            double *lv_J_kg)
{
  if (!table || !p_sat_Pa || !lv_J_kg || !isfinite(temperature_K) ||
      temperature_K <= 0. || table->n < 2)
    return P1P_ERR_ARGUMENT;
  if (temperature_K < table->temperature_K[0] ||
      temperature_K > table->temperature_K[table->n - 1])
    return P1P_ERR_RANGE;
  size_t hi = 1;
  while (hi < table->n && temperature_K > table->temperature_K[hi])
    hi++;
  if (temperature_K == table->temperature_K[table->n - 1])
    hi = table->n - 1;
  const size_t lo = hi - 1;
  const double f = (temperature_K - table->temperature_K[lo]) /
                   (table->temperature_K[hi] - table->temperature_K[lo]);
  const double logp = log(table->saturation_pressure_Pa[lo]) + f *
      (log(table->saturation_pressure_Pa[hi]) -
       log(table->saturation_pressure_Pa[lo]));
  *p_sat_Pa = exp(logp);
  *lv_J_kg = table->vaporization_enthalpy_J_kg[lo] + f *
      (table->vaporization_enthalpy_J_kg[hi] -
       table->vaporization_enthalpy_J_kg[lo]);
  if (!isfinite(*p_sat_Pa) || *p_sat_Pa <= 0. ||
      !isfinite(*lv_J_kg) || *lv_J_kg <= 0.)
    return P1P_ERR_TABLE;
  return P1P_OK;
}

/* Hertz-Knudsen free-molecular evaporation flux; negative net flux is clipped
 * to zero because P1P does not model condensation. R_specific is J/(kg K).
 */
static inline P1PStatus p1p_hk_mass_flux(double alpha_ev,
                                         double p_sat_Pa,
                                         double p_ag_vapor_Pa,
                                         double temperature_K,
                                         double r_ag_specific_J_kg_K,
                                         double *flux_kg_m2_s)
{
  if (!flux_kg_m2_s || !isfinite(alpha_ev) || alpha_ev < 0. || alpha_ev > 1. ||
      !isfinite(p_sat_Pa) || p_sat_Pa < 0. ||
      !isfinite(p_ag_vapor_Pa) || p_ag_vapor_Pa < 0. ||
      !isfinite(temperature_K) || temperature_K <= 0. ||
      !isfinite(r_ag_specific_J_kg_K) || r_ag_specific_J_kg_K <= 0.)
    return P1P_ERR_ARGUMENT;
  const double delta_p = fmax(0., p_sat_Pa - p_ag_vapor_Pa);
  const double denominator =
      sqrt(2. * P1P_PI * r_ag_specific_J_kg_K * temperature_K);
  if (!isfinite(denominator) || denominator <= 0.)
    return P1P_ERR_ARGUMENT;
  *flux_kg_m2_s = alpha_ev * delta_p / denominator;
  if (!isfinite(*flux_kg_m2_s) || *flux_kg_m2_s < 0.)
    return P1P_ERR_ARGUMENT;
  return P1P_OK;
}

/* Convert a surface flux into a local timestep update, capped by available
 * liquid Ag mass. The returned energy is latent-only: the caller must remove
 * the corresponding condensed-material enthalpy consistently when updating
 * composition, so the total escaped energy is not double-counted.
 */
static inline P1PStatus p1p_evaporation_step(double flux_kg_m2_s,
                                            double interface_area_m2,
                                            double dt_s,
                                            double available_liquid_mass_kg,
                                            double lv_J_kg,
                                            P1PStep *step)
{
  if (!step || !isfinite(flux_kg_m2_s) || flux_kg_m2_s < 0. ||
      !isfinite(interface_area_m2) || interface_area_m2 < 0. ||
      !isfinite(dt_s) || dt_s <= 0. ||
      !isfinite(available_liquid_mass_kg) || available_liquid_mass_kg < 0. ||
      !isfinite(lv_J_kg) || lv_J_kg < 0.)
    return P1P_ERR_ARGUMENT;
  step->mass_flux_kg_m2_s = flux_kg_m2_s;
  step->requested_mass_kg = flux_kg_m2_s * interface_area_m2 * dt_s;
  if (!isfinite(step->requested_mass_kg))
    return P1P_ERR_ARGUMENT;
  step->removed_mass_kg = fmin(step->requested_mass_kg,
                               available_liquid_mass_kg);
  step->latent_energy_J = step->removed_mass_kg * lv_J_kg;
  if (!isfinite(step->latent_energy_J))
    return P1P_ERR_ARGUMENT;
  step->mass_limited = step->requested_mass_kg > available_liquid_mass_kg;
  return P1P_OK;
}

/* Apply one cell-local, conservative surface source.  interface_area_m2 must
 * already be the PLIC area exposed to gas and weighted by the local Ag surface
 * coverage and liquid fraction.  The same removed Ag volume is removed from
 * both the composite VOF sample fraction and the Ag component fraction; the
 * fixed skeleton and pore fractions are unchanged.  The removed condensed Ag
 * sensible enthalpy and vaporization enthalpy are both booked as escaped
 * energy.  A local enthalpy floor of zero (300 K reference) prevents an
 * unphysical below-reference state; heat conducted/source-added in the next
 * thermal solve can replenish the cell in later steps. */
static inline P1PStatus p1p_apply_cell_evaporation(
    const P1PPropertyTable *table, double alpha_ev, double p_ag_vapor_Pa,
    double r_ag_specific_J_kg_K, double temperature_K,
    double liquid_fraction,
    double ag_liquid_specific_enthalpy_J_kg, double rho_ag_kg_m3,
    double cell_volume_m3, double interface_area_m2, double dt_s,
    P1PCellState *state, P1PCellUpdate *update)
{
  if (!table || !state || !update || !isfinite(state->sample_fraction) ||
      !isfinite(state->ag_fraction) || !isfinite(state->skeleton_fraction) ||
      !isfinite(state->pore_fraction) ||
      !isfinite(state->enthalpy_density_J_m3) ||
      !isfinite(liquid_fraction) || liquid_fraction < 0. ||
      liquid_fraction > 1. ||
      !isfinite(ag_liquid_specific_enthalpy_J_kg) ||
      ag_liquid_specific_enthalpy_J_kg < 0. || !isfinite(rho_ag_kg_m3) ||
      rho_ag_kg_m3 <= 0. || !isfinite(cell_volume_m3) || cell_volume_m3 <= 0. ||
      state->sample_fraction < 0. || state->sample_fraction > 1. ||
      state->ag_fraction < 0. || state->skeleton_fraction < 0. ||
      state->pore_fraction < 0. ||
      fabs(state->sample_fraction - state->ag_fraction -
           state->skeleton_fraction - state->pore_fraction) > 1.e-8 ||
      state->enthalpy_density_J_m3 < 0.)
    return P1P_ERR_ARGUMENT;

  memset(update, 0, sizeof(*update));
  double p_sat_Pa, lv_J_kg, flux_kg_m2_s;
  P1PStatus status = p1p_property_lookup(table, temperature_K,
                                         &p_sat_Pa, &lv_J_kg);
  if (status != P1P_OK)
    return status;
  double alpha_from_table = NAN;
  int used_table_alpha = p1p_alpha_lookup(table,temperature_K,
                                           &alpha_from_table) == P1P_OK;
  const double alpha_used = used_table_alpha ? alpha_from_table : alpha_ev;
  status = p1p_hk_mass_flux(alpha_used, p_sat_Pa, p_ag_vapor_Pa,
                            temperature_K, r_ag_specific_J_kg_K,
                            &flux_kg_m2_s);
  if (status != P1P_OK)
    return status;
  update->alpha_ev_used = alpha_used;
  update->alpha_from_property_table = used_table_alpha;

  const double available_liquid_mass_kg =
      state->ag_fraction * liquid_fraction * rho_ag_kg_m3 * cell_volume_m3;
  if (!isfinite(available_liquid_mass_kg))
    return P1P_ERR_ARGUMENT;
  /* The caller supplies an area already weighted by liquid fraction; this
   * matching cap limits the source to represented liquid Ag mass. */
  P1PStep step;
  status = p1p_evaporation_step(flux_kg_m2_s, interface_area_m2, dt_s,
                                available_liquid_mass_kg, lv_J_kg, &step);
  if (status != P1P_OK)
    return status;

  const double vapor_specific_energy_J_kg =
      ag_liquid_specific_enthalpy_J_kg + lv_J_kg;
  if (!isfinite(vapor_specific_energy_J_kg) ||
      !(vapor_specific_energy_J_kg > 0.))
    return P1P_ERR_ARGUMENT;
  const double energy_available_J =
      state->enthalpy_density_J_m3 * cell_volume_m3;
  if (!isfinite(energy_available_J))
    return P1P_ERR_ARGUMENT;
  const double energy_limited_mass_kg =
      energy_available_J / vapor_specific_energy_J_kg;
  if (!isfinite(energy_limited_mass_kg) || energy_limited_mass_kg < 0.)
    return P1P_ERR_ARGUMENT;
  if (step.removed_mass_kg > energy_limited_mass_kg) {
    step.removed_mass_kg = energy_limited_mass_kg;
    step.latent_energy_J = step.removed_mass_kg * lv_J_kg;
    update->enthalpy_limited = 1;
  }
  update->liquid_mass_limited = step.mass_limited;

  const double removed_volume_fraction =
      step.removed_mass_kg / (rho_ag_kg_m3 * cell_volume_m3);
  if (removed_volume_fraction > state->ag_fraction + 1.e-12 ||
      removed_volume_fraction > state->sample_fraction + 1.e-12)
    return P1P_ERR_ARGUMENT;
  const double sensible_energy_J =
      step.removed_mass_kg * ag_liquid_specific_enthalpy_J_kg;
  const double escaped_energy_J = sensible_energy_J + step.latent_energy_J;
  if (!isfinite(removed_volume_fraction) || removed_volume_fraction < 0. ||
      !isfinite(sensible_energy_J) || !isfinite(escaped_energy_J))
    return P1P_ERR_ARGUMENT;

  /* Validate a candidate update before mutating caller state.  A rejected
   * source must not leave a partially evaporated cell behind. */
  const double candidate_ag = state->ag_fraction - removed_volume_fraction;
  const double candidate_sample = state->sample_fraction - removed_volume_fraction;
  const double candidate_enthalpy =
      state->enthalpy_density_J_m3 - escaped_energy_J / cell_volume_m3;
  if (!isfinite(candidate_ag) || !isfinite(candidate_sample) ||
      !isfinite(candidate_enthalpy))
    return P1P_ERR_ARGUMENT;
  double ag_after = candidate_ag;
  double sample_after = candidate_sample;
  double enthalpy_after = candidate_enthalpy;
  /* Remove roundoff at a saturated mass/energy cap, never hide a physical
   * violation by clipping component fractions or enthalpy. */
  if (ag_after < 0. && ag_after > -1.e-12)
    ag_after = 0.;
  if (sample_after < 0. && sample_after > -1.e-12)
    sample_after = 0.;
  if (enthalpy_after < 0. && enthalpy_after > -1.e-7)
    enthalpy_after = 0.;
  if (ag_after < 0. || sample_after < 0. || enthalpy_after < 0. ||
      fabs(sample_after - ag_after -
           state->skeleton_fraction - state->pore_fraction) > 1.e-8)
    return P1P_ERR_ARGUMENT;

  state->ag_fraction = ag_after;
  state->sample_fraction = sample_after;
  state->enthalpy_density_J_m3 = enthalpy_after;

  update->mass_flux_kg_m2_s = flux_kg_m2_s;
  update->requested_mass_kg = step.requested_mass_kg;
  update->removed_mass_kg = step.removed_mass_kg;
  update->ag_sensible_energy_J = sensible_energy_J;
  update->latent_energy_J = step.latent_energy_J;
  update->escaped_energy_J = escaped_energy_J;
  update->removed_volume_fraction = removed_volume_fraction;
  return P1P_OK;
}

#endif
