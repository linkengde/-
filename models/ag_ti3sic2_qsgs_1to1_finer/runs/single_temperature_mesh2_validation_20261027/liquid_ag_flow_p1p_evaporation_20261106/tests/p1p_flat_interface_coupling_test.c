#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../p1p_evaporation.h"

static int check(double actual, double expected, double tolerance)
{
  double error=fabs(actual-expected);
  if (expected != 0.) error/=fabs(expected);
  if (!isfinite(actual) || error > tolerance) {
    fprintf(stderr,"FAIL actual=%.17g expected=%.17g tolerance=%.9g\n",
            actual,expected,tolerance);
    return 1;
  }
  return 0;
}

int main(int argc, char **argv)
{
  if (argc != 3) {
    fprintf(stderr,"usage: %s synthetic_property_table.csv output.csv\n",argv[0]);
    return 2;
  }
  P1PPropertyTable table;
  if (p1p_load_property_table(argv[1],&table)!=P1P_OK) {
    fprintf(stderr,"cannot load synthetic test property table\n");
    return 2;
  }
  double p_sat,lv,flux;
  const double T=1500.,alpha=0.5,p_v=0.,R_Ag=77.1;
  const double rho_Ag=10490.,h_Ag_liquid=4.e5,coverage=0.8;
  if (p1p_property_lookup(&table,T,&p_sat,&lv)!=P1P_OK ||
      p1p_hk_mass_flux(alpha,p_sat,p_v,T,R_Ag,&flux)!=P1P_OK) {
    fprintf(stderr,"synthetic property lookup failed\n");
    return 2;
  }
  FILE *out=fopen(argv[2],"w");
  if (!out) { perror(argv[2]); return 2; }
  fputs("grid_N,dt_s,steps,interface_area_m2,expected_mass_kg,actual_mass_kg,mass_relative_error,Ag_mass_residual,energy_residual,status\n",out);
  int failures=0;
  const int grids[]={8,16,32,64};
  const int substeps[]={1,2,4};
  const double total_time_s=1.e-6;
  for (size_t g=0;g<sizeof(grids)/sizeof(grids[0]);g++)
    for (size_t s=0;s<sizeof(substeps)/sizeof(substeps[0]);s++) {
      const int N=grids[g],steps=substeps[s];
      const double dx=1./N,cell_volume=dx*dx*dx,cell_area=dx*dx;
      const double dt=total_time_s/steps;
      double expected_area=1.,expected_mass=flux*coverage*expected_area*total_time_s;
      double initial_mass=0.,final_mass=0.,removed_mass=0.;
      double initial_energy=0.,final_energy=0.,escaped_energy=0.;
      P1PCellState face_cells[64*64];
      for (int face_cell=0;face_cell<N*N;face_cell++) {
        face_cells[face_cell]=(P1PCellState){0.5,0.4,0.1,0.,1.e9};
        initial_mass+=face_cells[face_cell].ag_fraction*rho_Ag*cell_volume;
        initial_energy+=face_cells[face_cell].enthalpy_density_J_m3*cell_volume;
      }
      for (int step_index=0;step_index<steps;step_index++)
        for (int face_cell=0;face_cell<N*N;face_cell++) {
          P1PCellUpdate update;
          P1PStatus status=p1p_apply_cell_evaporation(
              &table,alpha,p_v,R_Ag,T,1.,h_Ag_liquid,rho_Ag,
              cell_volume,coverage*cell_area,dt,&face_cells[face_cell],&update);
          if (status!=P1P_OK) {
            fprintf(stderr,"FAIL grid=%d step=%d cell=%d status=%d\n",
                    N,step_index,face_cell,(int)status);
            failures++;
            continue;
          }
          removed_mass+=update.removed_mass_kg;
          escaped_energy+=update.escaped_energy_J;
        }
      for (int face_cell=0;face_cell<N*N;face_cell++) {
        final_mass+=face_cells[face_cell].ag_fraction*rho_Ag*cell_volume;
        final_energy+=face_cells[face_cell].enthalpy_density_J_m3*cell_volume;
      }
      const double mass_error=fabs(removed_mass-expected_mass)/expected_mass;
      const double mass_residual=(initial_mass-final_mass-removed_mass)/initial_mass;
      const double energy_residual=(initial_energy-final_energy-escaped_energy)/initial_energy;
      int row_failed=0;
      row_failed+=check(removed_mass,expected_mass,2.e-12);
      row_failed+=check(mass_residual,0.,1.e-10);
      row_failed+=check(energy_residual,0.,1.e-10);
      failures+=row_failed;
      fprintf(out,"%d,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s\n",
              N,dt,steps,expected_area,expected_mass,removed_mass,mass_error,
              mass_residual,energy_residual,row_failed ? "FAIL" : "PASS");
    }
  fclose(out);
  if (failures) {
    fprintf(stderr,"FAIL P1P flat-interface coupling: %d failed checks\n",failures);
    return 1;
  }
  puts("PASS: synthetic flat-interface PLIC-area source, mesh and timestep coupling checks");
  return 0;
}
