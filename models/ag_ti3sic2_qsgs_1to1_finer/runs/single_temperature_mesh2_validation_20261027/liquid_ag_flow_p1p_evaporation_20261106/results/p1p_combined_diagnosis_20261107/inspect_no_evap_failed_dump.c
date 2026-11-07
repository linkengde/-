/* Offline reader for the archived synthetic no-evaporation failure dump.
 * It restores fields only; it does not call run() or advance the solution. */
#define main p1m_embedded_main
#include "../../tests/solver_smoke/original_solver_control/p1m_thermal_iteration_diagnostic.c"
#undef main

int main (int argc, char **argv)
{
  if (argc != 4) {
    fprintf(stderr, "usage: %s failed.dump property_dir cell_audit.csv\n", argv[0]);
    return 2;
  }
  read_sc_property_table(argv[2]);
  if (!restore(file = argv[1], list = all)) {
    fprintf(stderr, "restore failed: %s\n", argv[1]);
    return 3;
  }

  double minH=HUGE, maxH=-HUGE, minTinv=HUGE, maxTinv=-HUGE;
  double minT=HUGE, maxT=-HUGE, minGuess=HUGE, maxGuess=-HUGE;
  double maxPartition=0., maxDeltaT=0., maxDeltaH=0., maxInvResidual=0.;
  long cells=0, hotCells=0, tinyHotCells=0;
  double deltaX=0., deltaY=0., deltaZ=0.;
  double deltaF=0., deltaAg=0., deltaSc=0., deltaPore=0.;
  double deltaHadv=0., deltaHvol=0., deltaTfield=0., deltaTinv=0.;
  double invX=0., invY=0., invZ=0., invH=0., invT=0., invResidual=0.;
  FILE *fp=fopen(argv[3], "w");
  if (!fp) { perror(argv[3]); return 4; }
  fprintf(fp,"x_nm,y_nm,z_nm,f,ag,sc,pore,Hadv_J_m3,Hvol_J_m3,Tfield_K,T_from_Hadv_K,partition_error,H_inversion_residual_J_m3,delta_T_K\n");

  foreach (serial) {
    double f0=f[], a=ag_marker[], s=sc_marker[], p=pore_marker[];
    double H=Hadv[];
    double tinv=temperature_from_enthalpy(H,a,s);
    double cv,kk,hm;
    material_state(tinv,a,s,&cv,&kk,&hm);
    double partition=fabs(a+s+p-f0), hres=hm-H;
    double dT=fabs(tinv-Tfield[]), dH=fabs(Hvol[]-H);
    minH=min(minH,H); maxH=max(maxH,H);
    minTinv=min(minTinv,tinv); maxTinv=max(maxTinv,tinv);
    minT=min(minT,Tfield[]); maxT=max(maxT,Tfield[]);
    minGuess=min(minGuess,Tguess[]); maxGuess=max(maxGuess,Tguess[]);
    maxPartition=max(maxPartition,partition);
    if (fabs(hres) > maxInvResidual) {
      maxInvResidual=fabs(hres); invX=x; invY=y; invZ=z;
      invH=H; invT=tinv; invResidual=hres;
    }
    if (dT > maxDeltaT) {
      maxDeltaT=dT; deltaX=x; deltaY=y; deltaZ=z;
      deltaF=f0; deltaAg=a; deltaSc=s; deltaPore=p;
      deltaHadv=H; deltaHvol=Hvol[]; deltaTfield=Tfield[]; deltaTinv=tinv;
    }
    maxDeltaH=max(maxDeltaH,dH);
    cells++;
    if (tinv > 2000.) hotCells++;
    if (tinv > 2000. && f0 < 1.e-3) tinyHotCells++;
    if (!isfinite(tinv) || tinv > 2000. ||
        (f0 < 1.e-3 && H > 1.e6) || partition > 1.e-7 || dT > 10.)
      fprintf(fp,"%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
        x*1.e9,y*1.e9,z*1.e9,f0,a,s,p,H,Hvol[],Tfield[],tinv,
        partition,hres,dT);
  }
  fclose(fp);
  fprintf(stderr,
    "RESTORED cells=%ld Hadv_minmax_J_m3=%.12g,%.12g Tinv_minmax_K=%.12g,%.12g Tfield_minmax_K=%.12g,%.12g Tguess_minmax_K=%.12g,%.12g hot_cells=%ld tiny_hot_cells=%ld partition_error_max=%.12g max_abs_Tinv_minus_Tfield_K=%.12g at_nm=(%.12g,%.12g,%.12g) cell_fractions=(%.12g,%.12g,%.12g,%.12g) Hadv_Hvol_J_m3=(%.12g,%.12g) Tfield_Tinv_K=(%.12g,%.12g) max_abs_Hvol_minus_Hadv_J_m3=%.12g max_abs_inversion_residual_J_m3=%.12g at_nm=(%.12g,%.12g,%.12g) Hadv_J_m3=%.12g Tinv_K=%.12g signed_residual_J_m3=%.12g output=%s\n",
    cells,minH,maxH,minTinv,maxTinv,minT,maxT,minGuess,maxGuess,hotCells,
    tinyHotCells,maxPartition,maxDeltaT,deltaX*1.e9,deltaY*1.e9,deltaZ*1.e9,
    deltaF,deltaAg,deltaSc,deltaPore,deltaHadv,deltaHvol,deltaTfield,deltaTinv,
    maxDeltaH,maxInvResidual,invX*1.e9,invY*1.e9,invZ*1.e9,invH,invT,
    invResidual,argv[3]);
  return 0;
}
