/*
 * P1G exploratory 3-D fixed-grid Ag/Ag70 free-surface prototype.
 *
 * This is intentionally isolated from the production enthalpy solver.  It
 * combines a single-temperature enthalpy update, Gaussian composite heating,
 * fixed Ti3SiC2 Brinkman penalization, mushy/solid Ag damping, CLSVOF, 3-D
 * height-function capillary curvature, and projected Marangoni stress. All
 * physical outputs are exploratory; the
 * one-voxel-scale porous geometry and dynamic contact line are unresolved.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof-p1g.h"
#include "iforce.h"
#include "curvature.h"
#include "p1h_arc_pressure.h"

#define NC 64
#ifndef P1H_RECT_NZ
#define P1H_RECT_NZ 64
#endif
#define NCZ P1H_RECT_NZ
#define NCELLS (NC*NC*NCZ)
#ifndef P1H_TOP_OPEN
#define P1H_TOP_OPEN 0
#endif
#ifndef POOL_REPORT_INTERVAL_PS
#define POOL_REPORT_INTERVAL_PS 0.05
#endif
#ifndef POOL_FIELD_INTERVAL_PS
#define POOL_FIELD_INTERVAL_PS 0.05
#endif
#ifndef POOL_CHECKPOINT_INTERVAL_PS
#define POOL_CHECKPOINT_INTERVAL_PS 0.2
#endif

static uint8_t material_code[NCELLS]; /* 0 gas, 1 Ag, 2 Ti3SiC2 */
static double source_map[NCELLS];
static uint8_t pool_mask[NCELLS], pool_seen[NCELLS];
static double pool_weight[NCELLS];
static int pool_queue[NCELLS];
/* C-order 64^3 field buffers for compact, post-processable snapshots. */
/* C-order fixed-grid field buffers, written atomically and converted to HDF5. */
#define SNAPSHOT_NFIELDS 20
static float pool_snapshot_fields[SNAPSHOT_NFIELDS][NCELLS];
#define P1H_PRESSURE_NFIELDS 4
static float p1h_pressure_snapshot_fields[P1H_PRESSURE_NFIELDS][NCELLS];
#define SC_TABLE_CAPACITY 256
static int sc_property_count = 0;
static double sc_property_T[SC_TABLE_CAPACITY];
static double sc_property_Cp[SC_TABLE_CAPACITY];
static double sc_property_k[SC_TABLE_CAPACITY];
static double sc_property_h[SC_TABLE_CAPACITY];
static double sc_property_step_K = 0.;
static char sc_property_path[1024];

static int sigma_case = -1; /* -1 reducing atmosphere, +1 oxygen-rich reference */
static int marangoni_enabled = 1;
static double marangoni_multiplier = 1.; /* scales only tangential thermal Marangoni acceleration */
static P1HArcPressure arc_pressure_cfg;
static int arc_pressure_configured = 0;
static double arc_pressure_max_effective_Pa = 0.;
static double arc_pressure_force_z_N = 0.;
static FILE * arc_pressure_history = NULL;
static int resumed_from_checkpoint = 0;
static int mesh_gcd = 1;
/* Read by the isolated CLSVOF init-hook overlay so a restart preserves the
 * conservative VOF field instead of reconstructing it from the redistanced
 * level set. */
int p1g_restart_mode = 0;
static char restart_path[1024];
static double restored_time_ps = 0.;
static double solid_drag_s = 1.e11;
static double max_time_s = 5.e-12;
static const double heat_target_K = 2100.;
static const double heat_hold_duration_s = 3.e-12;
static double heat_dt_max_s = 2.e-14;
static int heat_on = 1;
static int heat_hold_started = 0;
static double heat_hold_start_s = NAN, heat_hold_end_s = NAN;
static uint32_t heat_hold_snapshot_mask = 0;
static double heat_source_ramp_s = 0.;
static double heat_source_ramp_start_s = INFINITY;
static int heat_source_ramping = 0;
static int solidified_steps = 0;
static int balanced_steps = 0;
static const int heat_balance_required_steps = 10;
static double energy_initial = 0., input_energy = 0., sink_energy = 0.;
static double ag_volume_initial = 0., sc_volume_initial = 0.;
static double applied_source_scale = 1.;
static double max_temperature_seen = 300.;
#ifndef P1I_SOURCE_POWER_W
#define P1I_SOURCE_POWER_W 0.008931899111848596
#endif
static double source_power_target = P1I_SOURCE_POWER_W;
static const double source_map_reference_power = 0.008931899111848596;
static double next_pool_report_ps = POOL_REPORT_INTERVAL_PS;
static double next_pool_field_snapshot_ps = POOL_FIELD_INTERVAL_PS;
static double next_pool_checkpoint_ps = POOL_CHECKPOINT_INTERVAL_PS;
static FILE * history = NULL;
static double initial_surface_height[NC*NC];
static scalar Tfield[], Told[], Tguess[], Hvol[], Hadv[], cpvol[], kcell[], mucell[];
static scalar ag_marker[], sc_marker[], liquid_fraction[];
static scalar source_q[], rhs_heat[], lambda_heat[], sigmaf[], drag_cell[];
static scalar kappa_field[], phi_normal;
static scalar arc_pressure_raw[], arc_pressure_mask[], arc_pressure_effective[];
static scalar total_normal_stress_field[];
static vector marangoni_cell[];
static face vector kappa_heat[], viscous[];

const double Lbox = 1.931887901542944e-8;
const double Tm = 1235.08;
const double rho_Ag = 10490.;
const double rho_Ag_liquid = 9346.;
const double rho_SC = 4520.;
const double rho_air = 1.2;
const double cp_air = 1005.;
const double k_air = 0.026;
const double cp_Ag_s_300 = 25.4/0.1078682;
const double cp_Ag_s_Tm = 31.8/0.1078682;
const double cp_Ag_l = 33.47/0.1078682;
const double L_Ag = 104729.66082682385;
const double Ag_mushy_width_K = 10.; /* numerical regularization, not a measured range */
const double Ag_Tsol = 1230.08;
const double Ag_Tliq = 1240.08;
const double SC_Tmax = 2300.;
const double Tref = 300.;
const double delta_x = Lbox/NC;
const double source_x = 95.8832e-10;
const double source_y = 96.5943950771472e-10;

Tfield[back] = dirichlet(300.);
Tfield[front] = neumann(0.);
u.n[back] = dirichlet(0.);
u.t[back] = dirichlet(0.);
u.r[back] = dirichlet(0.);
#if P1H_TOP_OPEN
u.n[front] = neumann(0.); /* pressure outlet for the gas cap; Ag remains inside unless it reaches the top */
u.t[front] = neumann(0.);
u.r[front] = neumann(0.);
p[front] = dirichlet(0.); /* ambient gauge-pressure reference */
#else
u.n[front] = dirichlet(0.); /* closed gas-cap boundary: no Ag can leave */
u.t[front] = neumann(0.);
u.r[front] = neumann(0.);
p[front] = neumann(0.);
#endif

static inline double clampd (double x, double lo, double hi)
{
  return x < lo ? lo : x > hi ? hi : x;
}

static int gcd_positive (int a, int b)
{
  while (b) {
    int r = a % b;
    a = b;
    b = r;
  }
  return a;
}

static double env_double (const char * name, double fallback)
{
  const char * value = getenv(name);
  if (!value || !*value) return fallback;
  char * end = NULL;
  double parsed = strtod(value, &end);
  if (end == value || *end || !isfinite(parsed)) {
    fprintf(stderr,"invalid numeric environment value %s=%s\n",name,value);
    exit(2);
  }
  return parsed;
}

static void read_arc_pressure_config (void)
{
  memset(&arc_pressure_cfg,0,sizeof(arc_pressure_cfg));
  arc_pressure_cfg.enabled = 0;
  arc_pressure_cfg.peak_Pa = 0.;
  arc_pressure_cfg.center_x_m = source_x;
  arc_pressure_cfg.center_y_m = source_y;
  arc_pressure_cfg.radius_m = 5.e-9;
  arc_pressure_cfg.start_s = 0.;
  arc_pressure_cfg.ramp_s = 0.;
  arc_pressure_cfg.end_s = INFINITY;
  arc_pressure_cfg.period_x_m = Lbox;
  arc_pressure_cfg.period_y_m = Lbox;
  arc_pressure_cfg.profile = P1H_PRESSURE_UNIFORM;

  const char * enabled = getenv("P1H_ARC_PRESSURE_ENABLED");
  if (enabled && !strcmp(enabled,"1")) arc_pressure_cfg.enabled = 1;
  else if (enabled && strcmp(enabled,"0")) {
    fprintf(stderr,"P1H_ARC_PRESSURE_ENABLED must be 0 or 1\n"); exit(2);
  }
  arc_pressure_cfg.peak_Pa = env_double("P1H_ARC_PRESSURE_PEAK_PA",0.);
  arc_pressure_cfg.center_x_m = env_double("P1H_ARC_PRESSURE_CENTER_X_NM",source_x*1.e9)*1.e-9;
  arc_pressure_cfg.center_y_m = env_double("P1H_ARC_PRESSURE_CENTER_Y_NM",source_y*1.e9)*1.e-9;
  arc_pressure_cfg.radius_m = env_double("P1H_ARC_PRESSURE_RADIUS_NM",5.)*1.e-9;
  arc_pressure_cfg.start_s = env_double("P1H_ARC_PRESSURE_START_PS",0.)*1.e-12;
  arc_pressure_cfg.ramp_s = env_double("P1H_ARC_PRESSURE_RAMP_PS",0.02)*1.e-12;
  double end_ps = env_double("P1H_ARC_PRESSURE_END_PS",-1.);
  arc_pressure_cfg.end_s = end_ps < 0. ? HUGE : end_ps*1.e-12;
  const char * profile = getenv("P1H_ARC_PRESSURE_PROFILE");
  if (profile && !strcmp(profile,"periodic_gaussian"))
    arc_pressure_cfg.profile = P1H_PRESSURE_PERIODIC_GAUSSIAN;
  else if (profile && strcmp(profile,"uniform")) {
    fprintf(stderr,"P1H_ARC_PRESSURE_PROFILE must be uniform or periodic_gaussian\n"); exit(2);
  }
  if (arc_pressure_cfg.enabled && !(arc_pressure_cfg.peak_Pa > 0.)) {
    fprintf(stderr,"enabled arc pressure requires P1H_ARC_PRESSURE_PEAK_PA > 0\n"); exit(2);
  }
  if (arc_pressure_cfg.ramp_s < 0. ||
      (isfinite(arc_pressure_cfg.end_s) && arc_pressure_cfg.end_s <= arc_pressure_cfg.start_s) ||
      (arc_pressure_cfg.profile == P1H_PRESSURE_PERIODIC_GAUSSIAN &&
       !(arc_pressure_cfg.radius_m > 0.))) {
    fprintf(stderr,"invalid P1H arc-pressure timing/profile configuration\n"); exit(2);
  }
  heat_source_ramp_s = env_double("P1H_SOURCE_RAMP_PS",arc_pressure_cfg.ramp_s*1.e12)*1.e-12;
  if (heat_source_ramp_s < 0. ||
      (arc_pressure_cfg.enabled && fabs(heat_source_ramp_s-arc_pressure_cfg.ramp_s) > 1.e-18)) {
    fprintf(stderr,"P1H source and pressure off-ramp durations must match for synchronous shutdown\n"); exit(2);
  }
  arc_pressure_configured = 1;
  fprintf(stderr,"ARC_PRESSURE_CONFIG enabled=%d peak_Pa=%.9g profile=%s center_nm=(%.9g,%.9g) radius_nm=%.9g start_ps=%.9g ramp_ps=%.9g end_ps=%.9g source_off_ramp_ps=%.9g periodic_xy=1\n",
          arc_pressure_cfg.enabled,arc_pressure_cfg.peak_Pa,
          arc_pressure_cfg.profile == P1H_PRESSURE_UNIFORM ? "uniform" : "periodic_gaussian",
          arc_pressure_cfg.center_x_m*1.e9,arc_pressure_cfg.center_y_m*1.e9,
          arc_pressure_cfg.radius_m*1.e9,arc_pressure_cfg.start_s*1.e12,
          arc_pressure_cfg.ramp_s*1.e12,
          end_ps < 0. ? -1. : arc_pressure_cfg.end_s*1.e12,
          heat_source_ramp_s*1.e12);
}

static void ramp_arc_pressure_down_at (double cutoff_s)
{
  if (!arc_pressure_configured || !arc_pressure_cfg.enabled)
    return;
  double ramp_end = cutoff_s + arc_pressure_cfg.ramp_s;
  if (arc_pressure_cfg.ramp_s <= 0.)
    ramp_end = cutoff_s;
  if (!isfinite(arc_pressure_cfg.end_s) || ramp_end < arc_pressure_cfg.end_s)
    arc_pressure_cfg.end_s = ramp_end;
}

static void stop_heat_and_pressure_at (double cutoff_s)
{
  heat_on = 0;
  heat_source_ramp_start_s = cutoff_s;
  heat_source_ramping = heat_source_ramp_s > 0.;
  ramp_arc_pressure_down_at(cutoff_s);
}

static double heat_source_scale_at (double time_s)
{
  return p1h_heat_source_scale(heat_on,heat_source_ramping,
                               heat_source_ramp_start_s,heat_source_ramp_s,time_s);
}

static void open_arc_pressure_history (int append)
{
  arc_pressure_history = fopen("p1h_arc_pressure_history.csv",append ? "a" : "w");
  if (!arc_pressure_history) { perror("p1h_arc_pressure_history.csv"); exit(2); }
  fseek(arc_pressure_history,0,SEEK_END);
  if (ftell(arc_pressure_history) == 0) {
    fprintf(arc_pressure_history,
            "time_ps,heat_source_scale,configured_peak_Pa,max_effective_interface_Pa,integrated_pressure_force_z_N,pressure_enabled,pressure_profile\n");
    fflush(arc_pressure_history);
  }
}

static inline int index3 (int i, int j, int k)
{
  return (i*NC + j)*NCZ + k;
}

static inline int input_index3 (int i, int j, int k)
{
  return (i*NC + j)*NC + k;
}

static void read_raw_inputs (const char * dir)
{
  char path[1024];
  const size_t input_cells=(size_t)NC*NC*NC;
  uint8_t * input_material=malloc(input_cells);
  double * input_source=malloc(input_cells*sizeof(double));
  if (!input_material || !input_source) {
    fprintf(stderr,"cannot allocate base P1G input buffers\n"); exit(2);
  }
  snprintf (path, sizeof(path), "%s/material_code.u8", dir);
  FILE * fp = fopen (path, "rb");
  if (!fp) { perror(path); exit(2); }
  if (fread (input_material, 1, input_cells, fp) != input_cells) {
    fprintf(stderr, "short read: %s\n", path); exit(2);
  }
  fclose(fp);
  snprintf (path, sizeof(path), "%s/source_power_density_W_m3.f64", dir);
  fp = fopen (path, "rb");
  if (!fp) { perror(path); exit(2); }
  if (fread (input_source, sizeof(double), input_cells, fp) != input_cells) {
    fprintf(stderr, "short read: %s\n", path); exit(2);
  }
  fclose(fp);
  for (int i=0; i<NC; i++)
    for (int j=0; j<NC; j++)
      for (int k=0; k<NCZ; k++) {
        int id=index3(i,j,k);
        if (k<NC) {
          int in=input_index3(i,j,k);
          material_code[id]=input_material[in];
          source_map[id]=input_source[in];
        }
        else {
          material_code[id]=0;
          source_map[id]=0.;
        }
      }
  free(input_material);
  free(input_source);
  double sum = 0.;
  long count_ag=0, count_sc=0, count_gas=0;
  for (int n = 0; n < NCELLS; n++) {
    if (material_code[n] > 2) {
      fprintf(stderr, "invalid material code at %d\n", n); exit(2);
    }
    if (material_code[n] == 0 && source_map[n] != 0.) {
      fprintf(stderr, "nonzero heat source in gas at %d\n", n); exit(2);
    }
    sum += source_map[n];
    count_ag += material_code[n] == 1;
    count_sc += material_code[n] == 2;
    count_gas += material_code[n] == 0;
  }
  sum *= cube(delta_x);
  if (!(sum > 0.) ||
      fabs(sum - source_map_reference_power) > 2.e-6*source_map_reference_power) {
    fprintf(stderr, "source-map template integral %.12g W != expected %.12g W\n",
            sum, source_map_reference_power);
    exit(2);
  }
  fprintf(stderr, "INPUT source_template=%.12g W normalized_target=%.12g W grid=%dx%dx%d dx=%.9g m codes_ag=%ld codes_sc=%ld codes_gas=%ld\n", sum,source_power_target,NC,NC,NCZ,delta_x,count_ag,count_sc,count_gas);
}

static void read_sc_property_table (const char * dir)
{
  snprintf(sc_property_path,sizeof(sc_property_path),
           "%s/Ti3SiC2_thermophysical_fit_300-2300K.csv",dir);
  FILE * fp=fopen(sc_property_path,"r");
  if (!fp) { perror(sc_property_path); exit(2); }
  char line[2048];
  if (!fgets(line,sizeof(line),fp)) {
    fprintf(stderr,"empty Ti3SiC2 property table: %s\n",sc_property_path);
    exit(2);
  }
  while (fgets(line,sizeof(line),fp)) {
    char * field[20];
    int n=0;
    char * token=strtok(line,",\r\n");
    while (token && n<20) { field[n++]=token; token=strtok(NULL,",\r\n"); }
    if (n<15 || sc_property_count>=SC_TABLE_CAPACITY) {
      fprintf(stderr,"invalid Ti3SiC2 CSV row/capacity at row %d\n",sc_property_count+2);
      exit(2);
    }
    int i=sc_property_count++;
    sc_property_T[i]=strtod(field[0],NULL);
    sc_property_Cp[i]=strtod(field[2],NULL);
    sc_property_k[i]=strtod(field[4],NULL);
    sc_property_h[i]=strtod(field[13],NULL);
    if (!isfinite(sc_property_T[i]) || !isfinite(sc_property_Cp[i]) ||
        !isfinite(sc_property_k[i]) || !isfinite(sc_property_h[i]) ||
        sc_property_Cp[i]<=0. || sc_property_k[i]<=0.) {
      fprintf(stderr,"invalid Ti3SiC2 property value at row %d\n",i+2);
      exit(2);
    }
    if (i==1) sc_property_step_K=sc_property_T[1]-sc_property_T[0];
    if (i>0 && (sc_property_T[i]<=sc_property_T[i-1] ||
        fabs((sc_property_T[i]-sc_property_T[i-1])-sc_property_step_K)>1.e-8)) {
      fprintf(stderr,"Ti3SiC2 property temperatures must be uniformly increasing\n");
      exit(2);
    }
  }
  fclose(fp);
  if (sc_property_count<2 || fabs(sc_property_T[0]-300.)>1.e-8 ||
      fabs(sc_property_T[sc_property_count-1]-SC_Tmax)>1.e-8 ||
      fabs(sc_property_h[0])>1.e-8) {
    fprintf(stderr,"Ti3SiC2 table must cover 300--%.0f K with h(300)=0\n",SC_Tmax);
    exit(2);
  }
  fprintf(stderr,"SC_PROPERTIES file=%s rows=%d range_K=%.1f..%.1f step_K=%.1f rho_used_kg_m3=%.1f heat_target_K=%.1f table_ceiling_K=%.1f\n",
          sc_property_path,sc_property_count,sc_property_T[0],
          sc_property_T[sc_property_count-1],sc_property_step_K,rho_SC,
          heat_target_K,SC_Tmax);
}

static double sc_table_value (double T, const double * values)
{
  if (sc_property_count<2 || !(sc_property_step_K>0.)) return NAN;
  double tc=clampd(T,sc_property_T[0],sc_property_T[sc_property_count-1]);
  double q=(tc-sc_property_T[0])/sc_property_step_K;
  int lo=(int)floor(q);
  if (lo>=sc_property_count-1) return values[sc_property_count-1];
  lo=max(lo,0);
  double f=q-lo;
  return values[lo]+f*(values[lo+1]-values[lo]);
}

static double ag_liquid_fraction (double T)
{
  return clampd ((T - Ag_Tsol)/(Ag_Tliq - Ag_Tsol), 0., 1.);
}

static double ag_cp_app (double T)
{
  double slope = (cp_Ag_s_Tm - cp_Ag_s_300)/(Tm - Tref);
  if (T < Ag_Tsol)
    return cp_Ag_s_300 + slope*(T - Tref);
  if (T <= Ag_Tliq)
    return cp_Ag_s_Tm + L_Ag/Ag_mushy_width_K;
  return cp_Ag_l;
}

static double ag_h (double T)
{
  double slope = (cp_Ag_s_Tm - cp_Ag_s_300)/(Tm - Tref);
  double t = max(T, Tref);
  double h_solidus = cp_Ag_s_300*(Ag_Tsol - Tref) +
    0.5*slope*sq(Ag_Tsol - Tref);
  if (t < Ag_Tsol)
    return cp_Ag_s_300*(t - Tref) + 0.5*slope*sq(t - Tref);
  if (t <= Ag_Tliq)
    return h_solidus + cp_Ag_s_Tm*(t - Ag_Tsol) +
      L_Ag*(t - Ag_Tsol)/Ag_mushy_width_K;
  return h_solidus + cp_Ag_s_Tm*Ag_mushy_width_K + L_Ag +
    cp_Ag_l*(t - Ag_Tliq);
}

static double sc_cp (double T)
{
  return sc_table_value(T,sc_property_Cp);
}

static double sc_k (double T)
{
  return sc_table_value(T,sc_property_k);
}

static double sc_h (double T)
{
  return sc_table_value(T,sc_property_h);
}

static double gamma_ag (double T)
{
  if (sigma_case < 0) {
    double tc = clampd(T, 1219., 1734.);
    return (984. - 0.298*(tc - 1234.))*1.e-3;
  }
  double tc = clampd(T, 1350., 1520.);
  return (698. + 0.546*(tc - 1234.))*1.e-3;
}

static double mu_ag (double T)
{
  double tc = clampd(T, 1273., 1673.);
  /* NIST ThermoML Arrhenius fit; experiment range 1273--1673 K. */
  return exp(-0.4631 + 2089./tc)*1.e-3;
}

static void material_state (double T, double agfrac, double scfrac,
                            double * cv, double * conductivity, double * hden)
{
  agfrac = clampd(agfrac, 0., 1.);
  scfrac = clampd(scfrac, 0., 1. - agfrac);
  double airfrac = max(0., 1. - agfrac - scfrac);
  *cv = agfrac*rho_Ag*ag_cp_app(T) + scfrac*rho_SC*sc_cp(T) + airfrac*rho_air*cp_air;
  *conductivity = agfrac*(T < Tm ? 429. + (175. - 429.)*(T-Tref)/(Tm-Tref) : 175.) +
    scfrac*sc_k(T) + airfrac*k_air;
  *hden = agfrac*rho_Ag*ag_h(T) + scfrac*rho_SC*sc_h(T) +
    airfrac*rho_air*cp_air*(T - Tref);
}

/* Invert the cell-mixture volumetric enthalpy after conservative advection. */
static double temperature_from_enthalpy (double target_h, double agfrac, double scfrac)
{
  if (!isfinite(target_h) || target_h <= 0.) return Tref;
  double lo = Tref, hi = 5000.;
  double cv, kk, hh;
  material_state(hi,agfrac,scfrac,&cv,&kk,&hh);
  if (target_h >= hh) return hi;
  /* 24 bisections bound the temperature inversion error below 3.0e-4 K. */
  for (int it = 0; it < 24; it++) {
    double mid = 0.5*(lo + hi);
    material_state(mid,agfrac,scfrac,&cv,&kk,&hh);
    if (hh < target_h) lo = mid;
    else hi = mid;
  }
  return 0.5*(lo + hi);
}

static double total_enthalpy (void)
{
  double E = 0.;
  foreach(reduction(+:E)) {
    double scfrac = clampd(sc_marker[], 0., 1.);
    double agfrac = clampd(ag_marker[], 0., 1. - scfrac);
    double cv, kk, hh;
    material_state(Tfield[], agfrac, scfrac, &cv, &kk, &hh);
    E += hh*dv();
  }
  return E;
}

/* Measure the largest face-connected, majority-liquid Ag component at the
 * saved peak field. x/y neighbors wrap to honor the periodic boundaries. */
static void write_peak_pool_stats (double time_ps, double Tmax_SC)
{
  memset(pool_mask,0,sizeof(pool_mask));
  memset(pool_seen,0,sizeof(pool_seen));
  foreach() {
    int ii=clamp((int)((x-X0)/Delta),0,NC-1);
    int jj=clamp((int)((y-Y0)/Delta),0,NC-1);
    int kk=clamp((int)((z-Z0)/Delta),0,NCZ-1);
    int id=index3(ii,jj,kk);
    double theta=clampd(f[],0.,1.)*clampd(ag_marker[],0.,1.)*
      ag_liquid_fraction(Tfield[]);
    theta=clampd(theta,0.,1.);
    pool_weight[id]=theta;
    pool_mask[id]=theta>=0.5;
  }
  double cell_nm3=pow(delta_x*1.e9,3.);
  double best_volume=0.;
  int best_cells=0, best_min_i=0, best_max_i=-1;
  int best_min_j=0, best_max_j=-1, best_min_k=0, best_max_k=-1;
  int component_count=0;
  for (int id=0; id<NCELLS; id++) {
    if (!pool_mask[id] || pool_seen[id]) continue;
    component_count++;
    int head=0, tail=0, cells=0;
    double volume=0.;
    int min_i=NC, max_i=-1, min_j=NC, max_j=-1, min_k=NCZ, max_k=-1;
    pool_queue[tail++]=id;
    pool_seen[id]=1;
    while (head<tail) {
      int cur=pool_queue[head++];
      int kk=cur%NCZ, jj=(cur/NCZ)%NC, ii=cur/(NC*NCZ);
      cells++;
      volume+=pool_weight[cur]*cell_nm3;
      min_i=min(min_i,ii); max_i=max(max_i,ii);
      min_j=min(min_j,jj); max_j=max(max_j,jj);
      min_k=min(min_k,kk); max_k=max(max_k,kk);
      int ni[6], nj[6], nk[6];
      ni[0]=(ii+NC-1)%NC; nj[0]=jj; nk[0]=kk;
      ni[1]=(ii+1)%NC;    nj[1]=jj; nk[1]=kk;
      ni[2]=ii; nj[2]=(jj+NC-1)%NC; nk[2]=kk;
      ni[3]=ii; nj[3]=(jj+1)%NC;    nk[3]=kk;
      int nn=4;
      if (kk>0) { ni[nn]=ii; nj[nn]=jj; nk[nn++]=kk-1; }
      if (kk<NCZ-1) { ni[nn]=ii; nj[nn]=jj; nk[nn++]=kk+1; }
      for (int q=0; q<nn; q++) {
        int next=index3(ni[q],nj[q],nk[q]);
        if (pool_mask[next] && !pool_seen[next]) {
          pool_seen[next]=1;
          pool_queue[tail++]=next;
        }
      }
    }
    if (volume>best_volume) {
      best_volume=volume; best_cells=cells;
      best_min_i=min_i; best_max_i=max_i;
      best_min_j=min_j; best_max_j=max_j;
      best_min_k=min_k; best_max_k=max_k;
    }
  }
  FILE *fp=fopen("p1h_peak_pool_dimensions.csv","a+");
  if (!fp) { perror("p1h_peak_pool_dimensions.csv"); return; }
  fseek(fp,0,SEEK_END);
  if (ftell(fp)==0)
    fprintf(fp,"time_ps,Tmax_skeleton_K,connectivity,liquid_Ag_mask_threshold,components,largest_cells,largest_liquid_volume_nm3,x_min_nm,x_max_nm,x_size_nm,y_min_nm,y_max_nm,y_size_nm,z_min_nm,z_max_nm,z_size_nm\n");
  if (best_cells>0) {
    double dxnm=delta_x*1.e9;
    fprintf(fp,"%.9g,%.9g,6-face-periodic-xy,0.5,%d,%d,%.12g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
            time_ps,Tmax_SC,component_count,best_cells,best_volume,
            best_min_i*dxnm,(best_max_i+1)*dxnm,(best_max_i-best_min_i+1)*dxnm,
            best_min_j*dxnm,(best_max_j+1)*dxnm,(best_max_j-best_min_j+1)*dxnm,
            best_min_k*dxnm,(best_max_k+1)*dxnm,(best_max_k-best_min_k+1)*dxnm);
  }
  else
    fprintf(fp,"%.9g,%.9g,6-face-periodic-xy,0.5,%d,0,0,nan,nan,0,nan,nan,0,nan,nan,0\n",
            time_ps,Tmax_SC,component_count);
  fclose(fp);
}

/* Save the complete P1G diagnostic field set; each raw frame is atomic. */
static void write_pool_field_snapshot (double time_ps, const char * stage)
{
  foreach() {
    int ii=clamp((int)((x-X0)/Delta),0,NC-1);
    int jj=clamp((int)((y-Y0)/Delta),0,NC-1);
    int kk=clamp((int)((z-Z0)/Delta),0,NCZ-1);
    int id=index3(ii,jj,kk);
    pool_snapshot_fields[0][id]=(float)f[];
    pool_snapshot_fields[1][id]=(float)d[];
    pool_snapshot_fields[2][id]=(float)ag_marker[];
    pool_snapshot_fields[3][id]=(float)sc_marker[];
    pool_snapshot_fields[4][id]=(float)liquid_fraction[];
    pool_snapshot_fields[5][id]=(float)Tfield[];
    pool_snapshot_fields[6][id]=(float)Hvol[];
    pool_snapshot_fields[7][id]=(float)u.x[];
    pool_snapshot_fields[8][id]=(float)u.y[];
    pool_snapshot_fields[9][id]=(float)u.z[];
    pool_snapshot_fields[10][id]=(float)p[];
    pool_snapshot_fields[11][id]=(float)kappa_field[];
    pool_snapshot_fields[12][id]=(float)gamma_ag(Tfield[]);
    pool_snapshot_fields[13][id]=(float)sigmaf[];
    pool_snapshot_fields[14][id]=(float)marangoni_cell.x[];
    pool_snapshot_fields[15][id]=(float)marangoni_cell.y[];
    pool_snapshot_fields[16][id]=(float)marangoni_cell.z[];
    pool_snapshot_fields[17][id]=(float)source_q[];
    pool_snapshot_fields[18][id]=(float)drag_cell[];
    pool_snapshot_fields[19][id]=(float)material_code[id];
    p1h_pressure_snapshot_fields[0][id]=(float)arc_pressure_raw[];
    p1h_pressure_snapshot_fields[1][id]=(float)arc_pressure_mask[];
    p1h_pressure_snapshot_fields[2][id]=(float)arc_pressure_effective[];
    p1h_pressure_snapshot_fields[3][id]=(float)total_normal_stress_field[];
  }
  char final_path[160], temporary_path[176];
  snprintf(final_path,sizeof(final_path),"p1h_%s_%s_t%010.6f_ps.bin",
           marangoni_enabled ? "marangoni_on" : "marangoni_off",
           stage ? stage : "regular",time_ps);
  snprintf(temporary_path,sizeof(temporary_path),"%s.part",final_path);
  FILE *fp=fopen(temporary_path,"wb");
  if (!fp) { perror(temporary_path); return; }
  struct {
    char magic[8];
    uint32_t version, nx, ny, nz, nfields, case_flags;
    double time_ps, dx_m, x0_m, y0_m, z0_m;
  } header={{'P','1','H','S','N','A','P','1'},1,NC,NC,NCZ,SNAPSHOT_NFIELDS,
            (uint32_t)((marangoni_enabled ? 1 : 0) | (P1H_TOP_OPEN ? 2 : 0)),
            time_ps,delta_x,X0,Y0,Z0};
  int ok=fwrite(&header,sizeof(header),1,fp)==1 &&
    fwrite(pool_snapshot_fields,sizeof(pool_snapshot_fields),1,fp)==1;
  if (fclose(fp)!=0) ok=0;
  if (!ok || rename(temporary_path,final_path)!=0) {
    perror("write/rename pool field snapshot");
    remove(temporary_path);
    return;
  }
  fprintf(stderr,"P1H_FIELD_SNAPSHOT t_ps=%.9g stage=%s marangoni=%d fields=%d file=%s\n",
          time_ps,stage ? stage : "regular",marangoni_enabled,SNAPSHOT_NFIELDS,final_path);

  snprintf(final_path,sizeof(final_path),"p1h_pressure_%s_%s_t%010.6f_ps.bin",
           marangoni_enabled ? "marangoni_on" : "marangoni_off",
           stage ? stage : "regular",time_ps);
  snprintf(temporary_path,sizeof(temporary_path),"%s.part",final_path);
  fp=fopen(temporary_path,"wb");
  if (!fp) { perror(temporary_path); return; }
  struct {
    char magic[8];
    uint32_t version, nx, ny, nz, nfields, enabled;
    double time_ps, dx_m, x0_m, y0_m, z0_m;
  } pheader={{'P','1','H','P','R','S','0','1'},1,NC,NC,NCZ,P1H_PRESSURE_NFIELDS,
             (uint32_t)arc_pressure_cfg.enabled,time_ps,delta_x,X0,Y0,Z0};
  ok=fwrite(&pheader,sizeof(pheader),1,fp)==1 &&
    fwrite(p1h_pressure_snapshot_fields,sizeof(p1h_pressure_snapshot_fields),1,fp)==1;
  if (fclose(fp)!=0) ok=0;
  if (!ok || rename(temporary_path,final_path)!=0) {
    perror("write/rename P1H pressure field snapshot");
    remove(temporary_path);
    return;
  }
  fprintf(stderr,"P1H_PRESSURE_FIELD_SNAPSHOT t_ps=%.9g stage=%s fields=%d file=%s\n",
          time_ps,stage ? stage : "regular",P1H_PRESSURE_NFIELDS,final_path);
}

struct RestartState {
  char magic[8];
  uint32_t version;
  int32_t heat_on, solidified_steps, balanced_steps, sigma_case, case_flags;
  uint32_t nx, ny, nz;
  double energy_initial, input_energy, sink_energy;
  double ag_volume_initial, sc_volume_initial, applied_source_scale;
  double max_temperature_seen, next_pool_report_ps, next_pool_field_snapshot_ps;
  double time_ps;
  int32_t arc_pressure_enabled, arc_pressure_profile, heat_source_ramping;
  double arc_pressure_peak_Pa, arc_pressure_center_x_m, arc_pressure_center_y_m;
  double arc_pressure_radius_m, arc_pressure_start_s, arc_pressure_ramp_s;
  double arc_pressure_end_s;
  double heat_source_ramp_s, heat_source_ramp_start_s;
  double marangoni_multiplier, heat_hold_start_s, heat_hold_end_s;
  uint32_t heat_hold_started, heat_hold_snapshot_mask;
};

static void write_restart_checkpoint (double time_ps, int step)
{
  char final_dump[192], temp_dump[200], final_state[208], temp_state[216];
  char marker[224], temp_marker[232], pointer_tmp[96];
  snprintf(final_dump,sizeof(final_dump),"p1h_restart_i%010d_t%010.6f_ps.dump",step,time_ps);
  snprintf(temp_dump,sizeof(temp_dump),"%s.part",final_dump);
  snprintf(final_state,sizeof(final_state),"%s.state",final_dump);
  snprintf(temp_state,sizeof(temp_state),"%s.part",final_state);
  snprintf(marker,sizeof(marker),"%s.complete",final_dump);
  snprintf(temp_marker,sizeof(temp_marker),"%s.part",marker);
  dump(file=temp_dump);
  if (rename(temp_dump,final_dump)!=0) { perror("rename restart dump"); return; }
  struct RestartState state={{'P','1','I','R','S','T','1','\0'},1,
    heat_on,solidified_steps,balanced_steps,sigma_case,
    (marangoni_enabled ? 1 : 0) | (P1H_TOP_OPEN ? 2 : 0),
    NC,NC,NCZ,
    energy_initial,input_energy,sink_energy,ag_volume_initial,sc_volume_initial,
    applied_source_scale,max_temperature_seen,next_pool_report_ps,
    next_pool_field_snapshot_ps,time_ps,
    arc_pressure_cfg.enabled,(int32_t)arc_pressure_cfg.profile,heat_source_ramping,
    arc_pressure_cfg.peak_Pa,arc_pressure_cfg.center_x_m,arc_pressure_cfg.center_y_m,
    arc_pressure_cfg.radius_m,arc_pressure_cfg.start_s,arc_pressure_cfg.ramp_s,
    arc_pressure_cfg.end_s,heat_source_ramp_s,heat_source_ramp_start_s,
    marangoni_multiplier,heat_hold_start_s,heat_hold_end_s,
    (uint32_t)heat_hold_started,heat_hold_snapshot_mask};
  FILE *fp=fopen(temp_state,"wb");
  if (!fp) { perror(temp_state); return; }
  int ok=fwrite(&state,sizeof(state),1,fp)==1 &&
    fwrite(initial_surface_height,sizeof(initial_surface_height),1,fp)==1;
  if (fclose(fp)!=0) ok=0;
  if (!ok || rename(temp_state,final_state)!=0) {
    perror("write restart state"); remove(temp_state); return;
  }
  fp=fopen(temp_marker,"w");
  if (!fp) { perror(temp_marker); return; }
  fprintf(fp,"time_ps=%.12g\nstep=%d\n",time_ps,step);
  if (fclose(fp)!=0 || rename(temp_marker,marker)!=0) {
    perror("write restart completion marker"); remove(temp_marker); return;
  }
  snprintf(pointer_tmp,sizeof(pointer_tmp),"latest_checkpoint.txt.part");
  fp=fopen(pointer_tmp,"w");
  if (fp) {
    fprintf(fp,"%s\n",final_dump);
    fclose(fp);
    rename(pointer_tmp,"latest_checkpoint.txt");
  }
  fprintf(stderr,"RESTART_CHECKPOINT t_ps=%.9g step=%d file=%s\n",time_ps,step,final_dump);
}

static int read_restart_state (const char * dump_path)
{
  char state_path[1200], marker_path[1200];
  snprintf(state_path,sizeof(state_path),"%s.state",dump_path);
  snprintf(marker_path,sizeof(marker_path),"%s.complete",dump_path);
  FILE * fp=fopen(marker_path,"r");
  if (!fp) { fprintf(stderr,"missing complete restart marker: %s\n",marker_path); return 0; }
  fclose(fp);
  fp=fopen(state_path,"rb");
  if (!fp) { perror(state_path); return 0; }
  struct RestartState state;
  int ok=fread(&state,sizeof(state),1,fp)==1 &&
    fread(initial_surface_height,sizeof(initial_surface_height),1,fp)==1;
  fclose(fp);
  if (!ok || memcmp(state.magic,"P1IRST1",7) || state.version!=1) {
    fprintf(stderr,"invalid P1H restart sidecar: %s\n",state_path); return 0;
  }
  if (state.nx!=NC || state.ny!=NC || state.nz!=NCZ) {
    fprintf(stderr,"restart grid does not match requested rectangular dimensions\n"); return 0;
  }
  heat_on=state.heat_on;
  solidified_steps=state.solidified_steps;
  balanced_steps=state.balanced_steps;
  heat_hold_started=(int)state.heat_hold_started;
  heat_hold_start_s=state.heat_hold_start_s;
  heat_hold_end_s=state.heat_hold_end_s;
  heat_hold_snapshot_mask=state.heat_hold_snapshot_mask;
  energy_initial=state.energy_initial;
  input_energy=state.input_energy;
  sink_energy=state.sink_energy;
  ag_volume_initial=state.ag_volume_initial;
  sc_volume_initial=state.sc_volume_initial;
  applied_source_scale=state.applied_source_scale;
  max_temperature_seen=state.max_temperature_seen;
  next_pool_report_ps=state.next_pool_report_ps;
  next_pool_field_snapshot_ps=state.next_pool_field_snapshot_ps;
  restored_time_ps=state.time_ps;
  next_pool_checkpoint_ps = (floor((restored_time_ps + 1.e-9)/POOL_CHECKPOINT_INTERVAL_PS) + 1.)*
    POOL_CHECKPOINT_INTERVAL_PS;
  int expected_case_flags=(marangoni_enabled ? 1 : 0) | (P1H_TOP_OPEN ? 2 : 0);
  if (state.sigma_case!=sigma_case || state.case_flags!=expected_case_flags) {
    fprintf(stderr,"restart case does not match requested sigma/Marangoni/top-boundary mode\n"); return 0;
  }
  if (state.arc_pressure_enabled!=arc_pressure_cfg.enabled ||
      state.arc_pressure_profile!=(int32_t)arc_pressure_cfg.profile ||
      state.arc_pressure_peak_Pa!=arc_pressure_cfg.peak_Pa ||
      state.arc_pressure_center_x_m!=arc_pressure_cfg.center_x_m ||
      state.arc_pressure_center_y_m!=arc_pressure_cfg.center_y_m ||
      state.arc_pressure_radius_m!=arc_pressure_cfg.radius_m ||
      state.arc_pressure_start_s!=arc_pressure_cfg.start_s ||
      state.arc_pressure_ramp_s!=arc_pressure_cfg.ramp_s ||
      state.heat_source_ramp_s!=heat_source_ramp_s ||
      state.marangoni_multiplier!=marangoni_multiplier) {
    fprintf(stderr,"restart case does not match requested P1H pressure configuration\n"); return 0;
  }
  arc_pressure_cfg.end_s=state.arc_pressure_end_s;
  heat_source_ramping=state.heat_source_ramping;
  heat_source_ramp_start_s=state.heat_source_ramp_start_s;
  return 1;
}

static void refresh_properties (void)
{
  foreach() {
    double scfrac = clampd(sc_marker[], 0., 1.);
    ag_marker[] = clampd(ag_marker[], 0., 1. - scfrac);
    double agfrac = ag_marker[];
    liquid_fraction[] = agfrac > 0. ? ag_liquid_fraction(Tfield[]) : 0.;
    material_state(Tfield[], agfrac, scfrac, &cpvol[], &kcell[], &Hvol[]);
    double fl = liquid_fraction[];
    double drag = scfrac*solid_drag_s +
      agfrac*solid_drag_s*sq(1. - fl)/(cube(fl) + 1.e-3);
    drag_cell[] = drag;
    sigmaf[] = agfrac*liquid_fraction[]*gamma_ag(Tfield[]);
    mucell[] = agfrac > 1.e-8 ? mu_ag(Tfield[]) : 1.8e-5;
  }
  boundary ({Tfield,Hvol,cpvol,kcell,mucell,ag_marker,sc_marker,
             liquid_fraction,drag_cell,sigmaf,arc_pressure_raw,
             arc_pressure_mask,arc_pressure_effective,total_normal_stress_field});
  foreach_face() {
    double kl = max(kcell[], 1.e-8), kr = max(kcell[-1], 1.e-8);
    kappa_heat.x[] = 2.*kl*kr/(kl + kr);
    double ml = max(mucell[], 1.e-8), mr = max(mucell[-1], 1.e-8);
    viscous.x[] = max(1.e-8, 2.*ml*mr/(ml + mr));
  }
}

int main (int argc, char * argv[])
{
  if (argc < 3 || argc > 5) {
    fprintf(stderr,"usage: %s input_dir negative|positive [solid_drag_s^-1] [max_time_ps]\nSet P1H_MARANGONI=on|off, optional P1H_RESTART_DUMP=restart.dump, and P1H_ARC_PRESSURE_* values in the environment.\n",argv[0]);
    return 2;
  }
  read_raw_inputs(argv[1]);
  read_sc_property_table(argv[1]);
  if (!strcmp(argv[2],"negative")) sigma_case = -1;
  else if (!strcmp(argv[2],"positive")) sigma_case = 1;
  else return 2;
  if (argc >= 4) solid_drag_s = atof(argv[3]);
  if (argc >= 5) max_time_s = atof(argv[4])*1.e-12;
  const char * marangoni_env=getenv("P1H_MARANGONI");
  if (!marangoni_env) marangoni_env=getenv("P1G_MARANGONI");
  if (marangoni_env && !strcmp(marangoni_env,"off")) marangoni_enabled=0;
  else if (marangoni_env && strcmp(marangoni_env,"on")) {
    fprintf(stderr,"P1H_MARANGONI must be on or off\n"); return 2;
  }
  marangoni_multiplier=env_double("P1H_MARANGONI_MULTIPLIER",1.);
  if (marangoni_multiplier < 0.) {
    fprintf(stderr,"P1H_MARANGONI_MULTIPLIER must be non-negative\n"); return 2;
  }
  const char * restart_env=getenv("P1H_RESTART_DUMP");
  if (!restart_env) restart_env=getenv("P1G_RESTART_DUMP");
  if (restart_env && *restart_env) snprintf(restart_path,sizeof(restart_path),"%s",restart_env);
  if (solid_drag_s <= 0. || max_time_s <= 0.) return 2;
  read_arc_pressure_config();

  /* Factor the finest cell counts so multigrid retains any exact common
   * coarsening levels on a rectangular domain (e.g. 64x64x76 -> 16x16x19). */
  mesh_gcd = gcd_positive(NC,NCZ);
  dimensions (nx = NC/mesh_gcd, ny = NC/mesh_gcd, nz = NCZ/mesh_gcd);
  size (Lbox);
  origin (0.,0.,0.);
  N = NC;
  periodic (left);
  periodic (bottom);
  rho1 = rho_Ag_liquid;
  rho2 = rho_air;
  mu = viscous;
  /* Advect total volumetric enthalpy and Ag occupancy conservatively; derive
   * temperature from H after advection, then solve the conductive enthalpy step. */
  tracers = list_concat (tracers, {Hvol, ag_marker});
  TOLERANCE = 1.e-5;
  NITERMAX = 200; /* allow extra iterations for the high-contrast liquid/mushy solve */
  CFL = 0.2;
  DT = heat_dt_max_s = 2.e-14;
  run();
}

event init (t = 0)
{
  if (restart_path[0]) {
    if (!restore(file=restart_path)) {
      fprintf(stderr,"cannot restore checkpoint %s\n",restart_path); exit(2);
    }
    p1g_restart_mode = 1;
    if (!read_restart_state(restart_path)) exit(2);
    resumed_from_checkpoint=1;
    refresh_properties();
    open_arc_pressure_history(1);
    history=fopen("p1h_history.csv","a");
    if (!history) { perror("p1h_history.csv append"); exit(2); }
    fseek(history,0,SEEK_END);
    if (ftell(history)==0) {
      fprintf(history,"time_ps,heat_on_after_step,Tmax_K,Tmax_Ag_K,Tmax_skeleton_K,mean_Ag_T_K,max_velocity_m_s,Ag_volume_nm3,Ag_volume_drift_pct,liquid_Ag_volume_nm3,skeleton_volume_nm3,composite_volume_nm3,mean_surface_depression_nm,source_center_depression_nm,source_power_W,enthalpy_J,input_energy_J,bottom_sink_energy_J,energy_residual_J,pressure_iterations,pressure_residual,thermal_iterations,thermal_residual,thermal_nonlinear_iterations,thermal_nonlinear_change_K,viscosity_iterations,viscosity_residual,dt_fs,heat_scale,input_minus_bottom_mW,balance_steps\n");
      fflush(history);
    }
    write_pool_field_snapshot(restored_time_ps,"restart_loaded");
    double restored_ag=0., restored_sc=0., restored_sample=0., restored_pin=0.;
    foreach(reduction(+:restored_ag) reduction(+:restored_sc) reduction(+:restored_sample)
            reduction(+:restored_pin)) {
      restored_ag+=ag_marker[]*dv();
      restored_sc+=sc_marker[]*dv();
      restored_sample+=f[]*dv();
      restored_pin+=source_q[]*f[]*dv();
    }
    fprintf(stderr,"RESTARTED t_ps=%.9g event_i=%d heat_on=%d Ag_nm3=%.12g skeleton_nm3=%.12g sample_nm3=%.12g source_power_mW=%.12g\n",
            restored_time_ps,i,heat_on,restored_ag/cube(1.e-9),restored_sc/cube(1.e-9),
            restored_sample/cube(1.e-9),restored_pin*1.e3);
    return 0;
  }
  double dx = Lbox/NC;
  foreach() {
    int i = clamp((int)((x-X0)/dx),0,NC-1);
    int j = clamp((int)((y-Y0)/dx),0,NC-1);
    int k = clamp((int)((z-Z0)/dx),0,NCZ-1);
    int code = material_code[index3(i,j,k)];
    d[] = 53.*delta_x - z; /* CLSVOF builds f from this signed-distance plane. */
    ag_marker[] = code == 1 ? 1. : 0.;
    Tfield[] = 300.;
    p[] = 0.;
    u.x[] = u.y[] = u.z[] = 0.;
    source_q[] = source_map[index3(i,j,k)];
  }
  boundary ({d,ag_marker,Tfield,p,u,source_q});
}

/* The CLSVOF init event builds f from d during i=0.  At i=1, initialize the
 * material tracers from that populated VOF field before starting heat flow. */
event initialize_material_fields (i = 1)
{
  if (resumed_from_checkpoint) return 0;
  double init_ag=0., init_sc=0., init_sample=0.;
  foreach() {
    int ii=clamp((int)((x-X0)/Delta),0,NC-1);
    int jj=clamp((int)((y-Y0)/Delta),0,NC-1);
    int kk=clamp((int)((z-Z0)/Delta),0,NCZ-1);
    int code=material_code[index3(ii,jj,kk)];
    ag_marker[] = code == 1 ? f[] : 0.;
    sc_marker[] = code == 2 ? f[] : 0.;
    arc_pressure_raw[] = arc_pressure_mask[] = arc_pressure_effective[] = 0.;
    total_normal_stress_field[] = 0.;
  }
  boundary ({ag_marker,sc_marker});
  double source_integral = 0.;
  foreach(reduction(+:source_integral))
    source_integral += source_q[]*f[]*dv();
  if (!(source_integral > 0.)) {
    fprintf(stderr,"ERROR: source integral is zero after VOF initialization\n");
    exit(2);
  }
  applied_source_scale = source_power_target/source_integral;
  foreach()
    source_q[] *= applied_source_scale;
  boundary ({source_q});
  refresh_properties();
  foreach(reduction(+:init_ag) reduction(+:init_sc) reduction(+:init_sample)) {
    init_ag+=ag_marker[]*dv();
    init_sc+=sc_marker[]*dv();
    init_sample+=f[]*dv();
  }
  ag_volume_initial = init_ag;
  sc_volume_initial = init_sc;
  for (int ij=0; ij<NC*NC; ij++)
    initial_surface_height[ij] = -HUGE;
  foreach() {
    if (f[] > 0.5) {
      int ii=clamp((int)((x-X0)/Delta),0,NC-1);
      int jj=clamp((int)((y-Y0)/Delta),0,NC-1);
      int ij=ii*NC+jj;
      initial_surface_height[ij] = max(initial_surface_height[ij],z + 0.5*Delta);
    }
  }
  fprintf(stderr,"INIT_FIELDS ag_m3=%.12g sc_m3=%.12g sample_m3=%.12g\n",init_ag,init_sc,init_sample);
  energy_initial = total_enthalpy();
  history = fopen("p1h_history.csv","w");
  if (!history) { perror("p1h_history.csv"); exit(2); }
  open_arc_pressure_history(0);
  fprintf(history,"time_ps,heat_on_after_step,Tmax_K,Tmax_Ag_K,Tmax_skeleton_K,mean_Ag_T_K,max_velocity_m_s,Ag_volume_nm3,Ag_volume_drift_pct,liquid_Ag_volume_nm3,skeleton_volume_nm3,composite_volume_nm3,mean_surface_depression_nm,source_center_depression_nm,source_power_W,enthalpy_J,input_energy_J,bottom_sink_energy_J,energy_residual_J,pressure_iterations,pressure_residual,thermal_iterations,thermal_residual,thermal_nonlinear_iterations,thermal_nonlinear_change_K,viscosity_iterations,viscosity_residual,dt_fs,heat_scale,input_minus_bottom_mW,balance_steps\n");
  fprintf(stderr,"INIT case=%s marangoni=%d marangoni_tangential_multiplier=%.9g top_boundary=%s grid=%dx%dx%d multigrid_factors=%dx%dx%d levels=%d drag_s-1=%.6g cells=%ld dx_m=%.9g dt_max_ps=%.6g source_power_W=%.12g source_scale=%.9g\n",
          sigma_case < 0 ? "negative_dsigma_dT" : "positive_dsigma_dT",marangoni_enabled,marangoni_multiplier,
          P1H_TOP_OPEN ? "open_ambient_pressure" : "closed_no_penetration",NC,NC,NCZ,
          NC/mesh_gcd,NC/mesh_gcd,NCZ/mesh_gcd,(int)grid->depth,solid_drag_s,
          grid->tn,delta_x,DT*1.e12,source_power_target,applied_source_scale);
  write_pool_field_snapshot(0.,"initial");
#ifndef P1I_PREFLIGHT
  write_restart_checkpoint(0.,i);
#endif
}

event thermal_step (i++,last)
{
  if (i == 0 || !(dt > 0.))
    return 0;
  /* Hvol has already passed through Basilisk's conservative tracer advection.
   * Recover the mixture temperature before applying conduction and the source. */
  foreach() {
    Hadv[] = Hvol[];
    double scfrac=clampd(sc_marker[],0.,1.);
    double agfrac=clampd(ag_marker[],0.,1.-scfrac);
    Tfield[]=temperature_from_enthalpy(Hadv[],agfrac,scfrac);
  }
  boundary ({Hadv,Tfield});
  int heat_was_on = heat_on;
  double max_sc_now = 300.;
  double heat_step_scale = 1.; /* Keep the configured source power fixed. */
  if (heat_on) {
    foreach(reduction(max:max_sc_now)) {
      double scfrac = clampd(sc_marker[],0.,1.);
      if (scfrac > 1.e-8)
        max_sc_now = max(max_sc_now,Tfield[]);
    }
    if (!heat_hold_started && max_sc_now >= heat_target_K) {
      heat_hold_started=1;
      heat_hold_start_s=t;
      heat_hold_end_s=t+heat_hold_duration_s;
      heat_hold_snapshot_mask=0;
      fprintf(stderr,"HEAT_HOLD_STARTED t_ps=%.9g Tmax_skeleton_K=%.9g hold_duration_ps=%.9g hold_end_ps=%.9g\n",
              t*1.e12,max_sc_now,heat_hold_duration_s*1.e12,heat_hold_end_s*1.e12);
      write_peak_pool_stats(t*1.e12,max_sc_now);
      write_pool_field_snapshot(t*1.e12,"hold_start_2100K");
      write_restart_checkpoint(t*1.e12,i);
    }
  }
  double heat_source_scale = heat_source_scale_at(t + 0.5*dt);

  double pin = 0., psink = 0.;
  foreach(reduction(+:pin))
    pin += heat_step_scale*heat_source_scale*source_q[]*f[]*dv();
  input_energy += dt*pin;

  /* Solve rho [h(Tnew)-h(Told)]/dt = div(k grad Tnew) + Q with a damped
   * Newton linearization of enthalpy.  This resolves the steep apparent Cp
   * across the Ag mushy interval without replacing latent heat by a clipped
   * temperature jump.  Scaling by dt/cv_ref conditions the Poisson solve. */
  const double cv_ref = 2.5e6; /* J m^-3 K^-1, conditioning scale only */
  foreach() {
    Told[] = Tfield[];
    Tguess[] = Tfield[];
  }
  boundary ({Told,Tguess});
  mgstats heat_mg;
  int thermal_nl_i = 0, thermal_nl_converged = 0;
  double thermal_nl_change = HUGE;
  for (thermal_nl_i = 0; thermal_nl_i < 80; thermal_nl_i++) {
    foreach() {
      Tguess[] = Tfield[];
      double scfrac = clampd(sc_marker[],0.,1.);
      double agfrac = clampd(ag_marker[],0.,1. - scfrac);
      double cvg, kval, hg;
      material_state(Tguess[],agfrac,scfrac,&cvg,&kval,&hg);
      double cv_eff = max(cvg, 1.);
      cpvol[] = cv_eff;
      kcell[] = kval;
      Hvol[] = hg;
      rhs_heat[] = ((hg - Hadv[] - cv_eff*Tguess[])/cv_ref -
                    heat_step_scale*heat_source_scale*source_q[]*f[]*dt/cv_ref);
      lambda_heat[] = -cv_eff/cv_ref;
    }
    boundary ({cpvol,kcell,Hvol,rhs_heat,lambda_heat,Tguess});
    foreach_face() {
      double kl=max(kcell[],1.e-8), kr=max(kcell[-1],1.e-8);
      kappa_heat.x[]=2.*kl*kr/(kl+kr)*dt/cv_ref;
    }
    heat_mg = poisson (Tfield, rhs_heat, kappa_heat, lambda_heat);
    boundary ({Tfield});
    thermal_nl_change = 0.;
    foreach(reduction(max:thermal_nl_change))
      thermal_nl_change = max(thermal_nl_change,fabs(Tfield[] - Tguess[]));
    if (thermal_nl_change < 1.e-2) {
      thermal_nl_converged = 1;
      break;
    }
    /* Mild damping keeps the Newton step stable when it crosses a material
     * property knot.  The converged fixed point is unchanged. */
    foreach()
      Tfield[] = 0.5*Tfield[] + 0.5*Tguess[];
    boundary ({Tfield});
  }
  if (!thermal_nl_converged) {
    fprintf(stderr,"STOP nonlinear enthalpy solve failed at t_ps=%.9g, change_K=%.9g after %d iterations\n",
            t*1.e12,thermal_nl_change,thermal_nl_i+1);
    dump(file="p1h_failed_nonlinear.dump");
    if (history) { fflush(history); fclose(history); history = NULL; }
    return 1;
  }
  boundary ({Tfield});
  refresh_properties();

  foreach(reduction(+:psink)) {
    if (fabs(z - (Z0 + 0.5*Delta)) < 0.25*Delta) {
      int ii=clamp((int)((x-X0)/Delta),0,NC-1);
      int jj=clamp((int)((y-Y0)/Delta),0,NC-1);
      int kk=0;
      int code=material_code[index3(ii,jj,kk)];
      if (code != 0)
        psink += 2.*kcell[]*max(0.,Tfield[]-300.)/Delta*sq(Delta);
    }
  }
  sink_energy += dt*psink;

  double Tmax=300., Tmax_Ag=300., Tmax_SC=300., agvol=0., liquidvol=0., scvol=0., samplevol=0., agtemp=0., agweight=0.;
  double umax=0.;
  double surface_dep_sum=0., surface_count=0.;
  double source_dep_sum=0., source_dep_count=0.;
  double heights[NC*NC];
  for (int ii=0; ii<NC*NC; ii++) heights[ii]=-HUGE;
  foreach(reduction(max:Tmax) reduction(max:Tmax_Ag) reduction(max:Tmax_SC) reduction(+:agvol) reduction(+:liquidvol)
          reduction(+:scvol) reduction(+:samplevol) reduction(+:agtemp) reduction(+:agweight)
          reduction(max:umax)) {
    int ii=clamp((int)((x-X0)/Delta),0,NC-1);
    int jj=clamp((int)((y-Y0)/Delta),0,NC-1);
    Tmax=max(Tmax,Tfield[]);
    umax=max(umax,norm(u));
    samplevol+=f[]*dv();
    double a=clampd(ag_marker[],0.,1. - clampd(sc_marker[],0.,1.));
    double s=clampd(sc_marker[],0.,1.);
    agvol+=a*dv();
    liquidvol+=a*ag_liquid_fraction(Tfield[])*dv();
    agtemp+=Tfield[]*a*dv();
    agweight+=a*dv();
    if (a > 1.e-8) Tmax_Ag=max(Tmax_Ag,Tfield[]);
    scvol+=s*dv();
    if (s > 1.e-8) Tmax_SC=max(Tmax_SC,Tfield[]);
    if (f[] > 0.5) {
      int ij=ii*NC+jj;
      heights[ij]=max(heights[ij],z+0.5*Delta);
    }
  }
  for (int ii=0; ii<NC; ii++)
    for (int jj=0; jj<NC; jj++) {
      int ij=ii*NC+jj;
      if (heights[ij] <= -HUGE/2. || initial_surface_height[ij] <= -HUGE/2.)
        continue;
      double depression=initial_surface_height[ij] - heights[ij];
      surface_dep_sum += depression;
      surface_count += 1.;
      double xc=(ii + 0.5)*delta_x;
      double yc=(jj + 0.5)*delta_x;
      if (sq(xc-source_x) + sq(yc-source_y) <= sq(5.e-9)) {
        source_dep_sum += depression;
        source_dep_count += 1.;
      }
    }
  double depression_nm=surface_count>0. ? surface_dep_sum/surface_count*1.e9 : NAN;
  double source_depression_nm=source_dep_count>0. ? source_dep_sum/source_dep_count*1.e9 : NAN;
  double input_minus_bottom_mW=(pin-psink)*1.e3;
  if (heat_on && !heat_hold_started && fabs(pin-psink) <= 0.01e-3)
    balanced_steps++;
  else
    balanced_steps=0;
  /* If the temperature crosses the target during this thermal update, begin
     the prescribed hold at this measured state. Do not stop on input/sink
     balance during the hold. */
  if (heat_on && !heat_hold_started && Tmax_SC >= heat_target_K) {
    heat_hold_started=1;
    heat_hold_start_s=t;
    heat_hold_end_s=t+heat_hold_duration_s;
    heat_hold_snapshot_mask=0;
    fprintf(stderr,"HEAT_HOLD_STARTED t_ps=%.9g Tmax_skeleton_K=%.9g hold_duration_ps=%.9g hold_end_ps=%.9g\n",
            t*1.e12,Tmax_SC,heat_hold_duration_s*1.e12,heat_hold_end_s*1.e12);
    write_peak_pool_stats(t*1.e12,Tmax_SC);
    write_pool_field_snapshot(t*1.e12,"hold_start_2100K");
    write_restart_checkpoint(t*1.e12,i);
  }
  /* Before the 2100 K hold begins, allow the source to stop when its net
     power into the sample has equilibrated with the bottom sink.  Requiring
     ten consecutive steps avoids stopping on a transient numerical crossing.
     Once the prescribed 3 ps hold has started, this balance criterion is
     disabled and only the property/divergence guards can interrupt it. */
  if (heat_on && !heat_hold_started &&
      balanced_steps >= heat_balance_required_steps) {
    stop_heat_and_pressure_at(t);
    fprintf(stderr,"HEAT_OFF_INPUT_SINK_BALANCE t_ps=%.9g Tmax_skeleton_K=%.9g input_minus_bottom_mW=%.9g consecutive_steps=%d hold_started=0\n",
            t*1.e12,Tmax_SC,input_minus_bottom_mW,balanced_steps);
  }
  if (heat_on && heat_hold_started) {
    const double offsets_ps[4]={0.5,1.,2.,3.};
    const char * stages[3]={"hold_plus_0p5ps","hold_plus_1ps","hold_plus_2ps"};
    for (int j=0; j<4; j++) {
      uint32_t bit=1u<<j;
      double event_s=heat_hold_start_s+offsets_ps[j]*1.e-12;
      if (!(heat_hold_snapshot_mask & bit) && t+1.e-18 >= event_s) {
        if (j < 3)
          write_pool_field_snapshot(t*1.e12,stages[j]);
        heat_hold_snapshot_mask |= bit;
        fprintf(stderr,"HEAT_HOLD_SNAPSHOT offset_ps=%.1f actual_t_ps=%.9g hold_age_ps=%.9g Tmax_skeleton_K=%.9g\n",
                offsets_ps[j],t*1.e12,(t-heat_hold_start_s)*1.e12,Tmax_SC);
        if (j < 3)
          write_restart_checkpoint(t*1.e12,i);
      }
    }
  }
  double rise_K = Tmax_SC - max_sc_now;
  double rate_K_s = rise_K > 0. && dt > 0. ? rise_K/dt : 0.;
  double ramp_peak_prediction_K = Tmax_SC +
    rate_K_s*heat_source_ramp_s*0.5;
  if (heat_on && Tmax_SC <= SC_Tmax &&
      rate_K_s > 0. && ramp_peak_prediction_K >= SC_Tmax &&
      (!heat_hold_started || t+1.e-18 < heat_hold_end_s)) {
    stop_heat_and_pressure_at(t);
    write_peak_pool_stats(t*1.e12,Tmax_SC);
    fprintf(stderr,"HEAT_OFF_PROPERTY_GUARD t_ps=%.9g Tmax_skeleton_K=%.9g predicted_after_ramp_K=%.9g Tmax_Ag_K=%.9g liquid_Ag_nm3=%.9g ceiling_K=%.9g hold_started=%d hold_age_ps=%.9g hold_complete=0\n",
            t*1.e12,Tmax_SC,ramp_peak_prediction_K,Tmax_Ag,
            liquidvol/cube(1.e-9),SC_Tmax,heat_hold_started,
            heat_hold_started ? (t-heat_hold_start_s)*1.e12 : 0.);
  }
  if (heat_on && heat_hold_started && t+1.e-18 >= heat_hold_end_s) {
    stop_heat_and_pressure_at(t);
    write_peak_pool_stats(t*1.e12,Tmax_SC);
    fprintf(stderr,"HEAT_OFF_AFTER_3PS_HOLD t_ps=%.9g Tmax_skeleton_K=%.9g Tmax_Ag_K=%.9g liquid_Ag_nm3=%.9g hold_age_ps=%.9g hold_complete=1\n",
            t*1.e12,Tmax_SC,Tmax_Ag,liquidvol/cube(1.e-9),
            (t-heat_hold_start_s)*1.e12);
  }
  /* Apply predictive timestep shrinking only while heating is still commanded
     on.  Once heat is switched off, the smooth source/pressure off-ramp is
     midpoint-integrated over the configured maximum step.  Keeping the
     temperature-rise cap active during this prescribed ramp caused a severe
     timestep collapse even while the measured peak temperature was falling.
     The hard 2300 K guard below remains active and never clips temperature. */
  double pending_heat_scale=heat_source_scale_at(t);
  if (heat_on && pending_heat_scale > 1.e-12 && Tmax_SC < SC_Tmax &&
      SC_Tmax - Tmax_SC < 20. && rate_K_s > 0.) {
    double safe_dt = 0.5*(SC_Tmax - Tmax_SC)/rate_K_s;
    DT = min(heat_dt_max_s,max(2.e-18,safe_dt));
  }
  else
    DT = heat_dt_max_s;
  if (heat_was_on && !heat_on) {
    write_pool_field_snapshot(t*1.e12,"heat_off");
    write_restart_checkpoint(t*1.e12,i);
  }
  double E=total_enthalpy();
  double residual=(E-energy_initial)-input_energy+sink_energy;
  double ag_drift_pct=ag_volume_initial>0. ? 100.*(agvol-ag_volume_initial)/ag_volume_initial : NAN;
  max_temperature_seen=max(max_temperature_seen,Tmax);
  if (history && (i%1==0)) {
  fprintf(history,"%.9g,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.12g,%.9g,%.12g,%.12g,%.12g,%.9g,%.9g,%.12g,%.12g,%.12g,%.12g,%.12g,%d,%.9g,%d,%.9g,%d,%.9g,%d,%.9g,%.9g,%.9g,%.9g,%d\n",
            t*1.e12,heat_on,Tmax,Tmax_Ag,Tmax_SC,agweight>0.?agtemp/agweight:300.,umax,
            agvol/cube(1.e-9),ag_drift_pct,liquidvol/cube(1.e-9),scvol/cube(1.e-9),samplevol/cube(1.e-9),depression_nm,
            source_depression_nm,
            pin,E,input_energy,sink_energy,residual,mgp.i,mgp.resa,heat_mg.i,heat_mg.resa,thermal_nl_i+1,
            thermal_nl_change,mgu.i,mgu.resa,dt*1.e15,heat_step_scale,input_minus_bottom_mW,balanced_steps);
    fflush(history);
  }
  if (i%10==0) {
    double ufmax=0.;
    foreach_face(reduction(max:ufmax))
      ufmax=max(ufmax,fabs(uf.x[]));
    fprintf(stderr,"STEP i=%d t_ps=%.6g Tmax=%.6gK Tmax_Ag=%.6gK Tmax_SC=%.6gK liquid_Ag_nm3=%.6g depression_nm=%.6g heat=%d dt_fs=%.9g dtmax_fs=%.9g ufmax_m_s=%.9g\n",
            i,t*1.e12,Tmax,Tmax_Ag,Tmax_SC,liquidvol/cube(1.e-9),depression_nm,heat_on,
            dt*1.e15,DT*1.e15,ufmax);
  }
  if (Tmax_SC>SC_Tmax) {
    fprintf(stderr,"STOP Ti3SiC2 property ceiling exceeded: %.9g K > %.9g K\n",Tmax_SC,SC_Tmax);
    dump(file="p1h_failed_property_ceiling.dump");
    if (history) { fflush(history); fclose(history); history = NULL; }
    return 1;
  }
  if (!isfinite(Tmax) || !isfinite(umax) || umax > 1000. || dt < 1.e-18) {
    fprintf(stderr,"STOP flow/thermal divergence guard: Tmax=%.9g K umax=%.9g m/s dt=%.9g s\n",
            Tmax,umax,dt);
    dump(file="p1h_failed_divergence.dump");
    if (history) { fflush(history); fclose(history); history = NULL; }
    return 1;
  }
  if (!heat_on && heat_source_scale_at(t) <= 1.e-12 &&
      liquidvol/cube(1.e-9) <= 1.e-8)
    solidified_steps++;
  else
    solidified_steps = 0;
  if (solidified_steps >= 2) {
    write_pool_field_snapshot(t*1.e12,"solidified");
    write_restart_checkpoint(t*1.e12,i);
    dump(file="p1h_final_solidified.dump");
    fprintf(stderr,"DONE_SOLIDIFIED t_ps=%.9g Tmax=%.9g K Tmax_Ag=%.9g K Tmax_skeleton=%.9g K\n",
            t*1.e12,Tmax,Tmax_Ag,Tmax_SC);
    return 1;
  }
}

/* Lightweight mid-run dimension reporting.  Use an iteration event and an
 * explicit physical-time threshold: Basilisk time events can be skipped when
 * dt adapts across their scheduled time in this coupled multiphysics run. */
event periodic_pool_dimensions (i++,last)
{
  double time_ps = t*1.e12;
  if (time_ps + 1.e-9 >= next_pool_report_ps) {
    double Tmax_SC=300.;
    foreach(reduction(max:Tmax_SC))
      if (sc_marker[]>1.e-8)
        Tmax_SC=max(Tmax_SC,Tfield[]);
    write_peak_pool_stats(time_ps,Tmax_SC);
    fprintf(stderr,"POOL_DIMENSIONS t_ps=%.9g Tmax_skeleton_K=%.9g\n",time_ps,Tmax_SC);
    if (time_ps + 1.e-9 >= next_pool_field_snapshot_ps) {
      write_pool_field_snapshot(time_ps,"regular");
      do {
        next_pool_field_snapshot_ps += POOL_FIELD_INTERVAL_PS;
      } while (next_pool_field_snapshot_ps <= time_ps + 1.e-9);
    }
    if (time_ps + 1.e-9 >= next_pool_checkpoint_ps) {
      write_restart_checkpoint(time_ps,i);
      do {
        next_pool_checkpoint_ps += POOL_CHECKPOINT_INTERVAL_PS;
      } while (next_pool_checkpoint_ps <= time_ps + 1.e-9);
    }
    do {
      next_pool_report_ps += POOL_REPORT_INTERVAL_PS;
    } while (next_pool_report_ps <= time_ps + 1.e-9);
  }
}

event acceleration (i++)
{
  /* Use Basilisk's 3-D height-function curvature for normal capillary
   * stress.  integral.h's distance-curvature/tensor implementation in this
   * pinned source is 2-D and must not be used on the octree grid. */
  cstats curvature_stats = curvature (f, kappa_field);
  (void) curvature_stats;
  phi_normal = new scalar;
  double pressure_max = 0.;
  foreach(reduction(max:pressure_max)) {
    double p_raw = p1h_arc_pressure_at(&arc_pressure_cfg,x,y,t);
    double p_mask = 0., p_effective = 0.;
    coord grad_f = {0.,0.,0.};
    marangoni_cell.x[] = marangoni_cell.y[] = marangoni_cell.z[] = 0.;
    if (f[] > 1.e-8 && f[] < 1. - 1.e-8) {
      grad_f = (coord){
        (f[1] - f[-1])/(2.*Delta),
        (f[0,1] - f[0,-1])/(2.*Delta),
        (f[0,0,1] - f[0,0,-1])/(2.*Delta)
      };
      coord grad_sigma = {
        (sigmaf[1] - sigmaf[-1])/(2.*Delta),
        (sigmaf[0,1] - sigmaf[0,-1])/(2.*Delta),
        (sigmaf[0,0,1] - sigmaf[0,0,-1])/(2.*Delta)
      };
      double gn = sqrt(sq(grad_f.x) + sq(grad_f.y) + sq(grad_f.z));
      if (gn > 1.e-12) {
        coord n = {grad_f.x/gn, grad_f.y/gn, grad_f.z/gn};
        double ag_surface_coverage = f[] > 1.e-12 ?
          clampd(ag_marker[]/f[],0.,1.) : 0.;
        p_mask = p1h_arc_pressure_mask(ag_surface_coverage,liquid_fraction[],n.z);
        p_effective = p_raw*p_mask;
        double normal_gradient = grad_sigma.x*n.x + grad_sigma.y*n.y + grad_sigma.z*n.z;
        double rho_mix = max(1., f[]*rho1 + (1. - f[])*rho2);
        marangoni_cell.x[] = (grad_sigma.x - normal_gradient*n.x)*gn/rho_mix;
        marangoni_cell.y[] = (grad_sigma.y - normal_gradient*n.y)*gn/rho_mix;
        marangoni_cell.z[] = (grad_sigma.z - normal_gradient*n.z)*gn/rho_mix;
        marangoni_cell.x[] *= marangoni_multiplier;
        marangoni_cell.y[] *= marangoni_multiplier;
        marangoni_cell.z[] *= marangoni_multiplier;
      }
    }
    double capillary_potential = kappa_field[] < nodata ?
      sigmaf[]*kappa_field[] : 0.;
    int interface_cell = f[] > 1.e-8 && f[] < 1. - 1.e-8;
    total_normal_stress_field[] = interface_cell ?
      capillary_potential + p_effective : 0.;
    phi_normal[] = interface_cell ? total_normal_stress_field[] : nodata;
    arc_pressure_raw[] = p_raw;
    arc_pressure_mask[] = p_mask;
    arc_pressure_effective[] = p_effective;
    pressure_max = max(pressure_max,p_effective);
  }
  boundary ({phi_normal,marangoni_cell,arc_pressure_raw,arc_pressure_mask,
             arc_pressure_effective,total_normal_stress_field});
  double pressure_force_z = 0.;
  foreach_face(z, reduction(+:pressure_force_z)) {
    int current_interface = f[] > 1.e-8 && f[] < 1. - 1.e-8;
    int lower_interface = f[0,0,-1] > 1.e-8 && f[0,0,-1] < 1. - 1.e-8;
    if (f[] != f[0,0,-1] && (current_interface || lower_interface)) {
      double p_face = current_interface && lower_interface ?
        0.5*(arc_pressure_effective[] + arc_pressure_effective[0,0,-1]) :
        current_interface ? arc_pressure_effective[] : arc_pressure_effective[0,0,-1];
      pressure_force_z += p_face*(f[] - f[0,0,-1])/Delta*dv();
    }
  }
  arc_pressure_force_z_N = pressure_force_z;
  arc_pressure_max_effective_Pa = pressure_max;
  f.phi = phi_normal;

  if (arc_pressure_history) {
    fprintf(arc_pressure_history,"%.12g,%.12g,%.12g,%.12g,%.12g,%d,%s\n",
            t*1.e12,heat_source_scale_at(t),arc_pressure_cfg.peak_Pa,arc_pressure_max_effective_Pa,
            arc_pressure_force_z_N,arc_pressure_cfg.enabled,
            arc_pressure_cfg.profile == P1H_PRESSURE_UNIFORM ? "uniform" : "periodic_gaussian");
    fflush(arc_pressure_history);
  }

  foreach_face() {
    double drag=0.5*(drag_cell[]+drag_cell[-1]);
    a.x[] += (marangoni_enabled ? 0.5*(marangoni_cell.x[] + marangoni_cell.x[-1]) : 0.) - drag*u.x[];
  }
}

event stop (t = max_time_s)
{
  double ufmax=0., umax=0.;
  foreach_face(reduction(max:ufmax))
    ufmax=max(ufmax,fabs(uf.x[]));
  foreach(reduction(max:umax))
    umax=max(umax,norm(u));
  fprintf(stderr,"STOP_DIAG t_ps=%.9g dt_fs=%.9g dtmax_fs=%.9g ufmax_m_s=%.9g umax_m_s=%.9g\n",
          t*1.e12,dt*1.e15,DT*1.e15,ufmax,umax);
  write_pool_field_snapshot(t*1.e12,"max_time");
#ifndef P1I_PREFLIGHT
  write_restart_checkpoint(t*1.e12,i);
  dump(file="p1h_final.dump");
#endif
  fprintf(stderr,"MAX_TIME_STOP t_ps=%.9g; check liquid_Ag_volume_nm3 in p1h_history.csv to determine whether cooling completed\n",
          t*1.e12);
  if (history) fclose(history);
  if (arc_pressure_history) { fclose(arc_pressure_history); arc_pressure_history = NULL; }
  return 1;
}
