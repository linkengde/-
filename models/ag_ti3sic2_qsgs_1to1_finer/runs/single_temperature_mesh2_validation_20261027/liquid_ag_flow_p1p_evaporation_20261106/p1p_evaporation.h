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
} P1PPropertyTable;

typedef struct {
  double mass_flux_kg_m2_s;
  double requested_mass_kg;
  double removed_mass_kg;
  double latent_energy_J;
  int mass_limited;
} P1PStep;

/* CSV contract: temperature_K,p_sat_Pa,L_v_J_kg; numeric rows only. */
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
  static const char *header = "temperature_K,p_sat_Pa,L_v_J_kg";
  if (strncmp(line, header, strlen(header)) != 0) {
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
    double t, p, lv;
    char tail;
    if (sscanf(line, " %lf , %lf , %lf %c", &t, &p, &lv, &tail) != 3 ||
        !isfinite(t) || !isfinite(p) || !isfinite(lv) ||
        t <= 0. || p <= 0. || lv <= 0. ||
        (table->n && t <= table->temperature_K[table->n - 1])) {
      fclose(fp);
      return P1P_ERR_TABLE;
    }
    size_t i = table->n++;
    table->temperature_K[i] = t;
    table->saturation_pressure_Pa[i] = p;
    table->vaporization_enthalpy_J_kg[i] = lv;
  }
  fclose(fp);
  return table->n >= 2 ? P1P_OK : P1P_ERR_TABLE;
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
  *flux_kg_m2_s = alpha_ev * delta_p /
      sqrt(2. * P1P_PI * r_ag_specific_J_kg_K * temperature_K);
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
  step->removed_mass_kg = fmin(step->requested_mass_kg,
                               available_liquid_mass_kg);
  step->latent_energy_J = step->removed_mass_kg * lv_J_kg;
  step->mass_limited = step->requested_mass_kg > available_liquid_mass_kg;
  return P1P_OK;
}

#endif
