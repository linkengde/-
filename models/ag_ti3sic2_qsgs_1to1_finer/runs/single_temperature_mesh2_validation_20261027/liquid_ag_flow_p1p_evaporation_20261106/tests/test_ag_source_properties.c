#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../p1p_evaporation.h"

static int failures = 0;
static FILE *out;

static void close_check(const char *name, double expected, double actual,
                        double relative_tolerance)
{
  double error = fabs(actual - expected)/fmax(fabs(expected), 1.e-300);
  int pass = isfinite(actual) && error <= relative_tolerance;
  fprintf(out, "%s,%.17g,%.17g,%.17g,%s\n", name, expected, actual,
          error, pass ? "PASS" : "FAIL");
  if (!pass) {
    fprintf(stderr,"FAIL %s expected=%.17g actual=%.17g relative_error=%.9g\n",
            name,expected,actual,error);
    failures++;
  }
}

static void status_check(const char *name, P1PStatus expected,
                         P1PStatus actual)
{
  int pass = expected == actual;
  fprintf(out,"%s,%d,%d,0,%s\n",name,(int)expected,(int)actual,
          pass ? "PASS" : "FAIL");
  if (!pass) failures++;
}

int main(int argc, char **argv)
{
  if (argc != 3) {
    fprintf(stderr,"usage: %s Ag_evaporation_solver_table.csv output.csv\n",argv[0]);
    return 2;
  }
  out=fopen(argv[2],"w");
  if (!out) { perror(argv[2]); return 2; }
  fputs("test,expected,actual,relative_error,status\n",out);
  P1PPropertyTable table;
  status_check("source_table_load",P1P_OK,p1p_load_property_table(argv[1],&table));
  double psat=NAN,lv=NAN,alpha=NAN;
  status_check("source_alpha_1373_K",P1P_OK,p1p_alpha_lookup(&table,1373.,&alpha));
  close_check("source_alpha_1373_value",0.2118,alpha,1.e-12);
  status_check("source_alpha_interpolation_1523_K",P1P_OK,
               p1p_alpha_lookup(&table,1523.,&alpha));
  close_check("source_alpha_1523_value",0.1131,alpha,1.e-12);
  status_check("source_alpha_upper_anchor_1773_K",P1P_OK,
               p1p_alpha_lookup(&table,1773.,&alpha));
  close_check("source_alpha_1773_value",0.2715,alpha,1.e-12);
  status_check("alpha_no_low_temperature_extrapolation",P1P_ERR_RANGE,
               p1p_alpha_lookup(&table,1372.99,&alpha));
  status_check("alpha_no_high_temperature_extrapolation",P1P_ERR_RANGE,
               p1p_alpha_lookup(&table,1773.01,&alpha));
  status_check("Wang_pressure_and_derived_Lv_2300_K",P1P_OK,
               p1p_property_lookup(&table,2300.,&psat,&lv));
  close_check("Wang_pressure_2300_Pa",50893.4361773,psat,1.e-10);
  close_check("derived_Lv_2300_J_kg",2405062.25,lv,1.e-8);
  status_check("table_alpha_precedes_sensitivity_fallback",P1P_OK,
               p1p_alpha_lookup(&table,1373.,&alpha));

  P1PCellState state={0.8,0.5,0.2,0.1,1.e8};
  P1PCellUpdate update;
  P1PStatus st=p1p_apply_cell_evaporation(&table,0.1,0.,77.1,1373.,1.,
      1.e6,10490.,1.e-12,1.e-12,1.e-12,&state,&update);
  status_check("real_source_table_cell_update",P1P_OK,st);
  close_check("cell_uses_table_alpha_not_fallback",0.2118,
              update.alpha_ev_used,1.e-12);
  close_check("cell_marks_experimental_alpha_source",1.,
              update.alpha_from_property_table,0.);
  fclose(out);
  if (failures) return 1;
  fprintf(stderr,"Ag source property checks passed: 12\n");
  return 0;
}
