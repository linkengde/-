# Reproduce the project's automatic-curvature P1b comparison

Set `BASILISK_SRC` to the `src` directory of the same pinned Basilisk snapshot recorded in `../p1b_official/p1b_source_manifest.json`, build `qcc` and the Tiny GL libraries, then run:

```sh
export BASILISK="$BASILISK_SRC"
export OPENGLIBS=-lfb_tiny
"$BASILISK_SRC/qcc" -disable-dimensions -O2 -Wall -autolink \
  p1b_vof_autocurvature.c -o p1b_vof_autocurvature -lm
timeout 900 ./p1b_vof_autocurvature \
  > p1b_vof_fullstress_stdout.log 2> p1b_vof_fullstress_stderr.log
```

`-disable-dimensions` is used because this is a fully nondimensional benchmark and the Basilisk dimensional checker otherwise reports its generic `timestep.h` constraints for the custom dimensionless fields. The implementation uses standard VOF (`two-phase.h`) and automatic curvature (`tension.h`), with a separately added CSF correction for spatially varying normal surface tension and its tangential gradient. It is a project-method comparison, not the official CLSVOF/integral-force discretization.
