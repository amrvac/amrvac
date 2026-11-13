#!/usr/bin/env bash
# Convergence study for phys='rho' (see the header of mod_usr.fpp).
#
# Builds once (every setting varied here is a runtime parameter), then runs
# every combination below and calls convergence.py for the table and plot.
# All generated par files, logs and stdout go to $OUTDIR; nothing is written
# in this directory except the build.
#
# Usage (from this directory, with AGILE_DIR set and the venv active):
#   ./run_convergence.sh                 # x sweep + diag sweep up to 200^3
#   WITH400=1 ./run_convergence.sh       # also diag at 400^3 (~13 GB, ~1 h on 4 cores)
#   SETUPS=x ./run_convergence.sh        # only the cheap quasi-1D sweep
# Other knobs (environment): NP (MPI ranks for diag, default 4), ARCH (gnu),
# OUTDIR, X_LIMITERS, X_STEPPERS, X_FLUXES, X_NLIST, D_LIMITERS, D_STEPPERS,
# D_FLUXES, D_NLIST, MPIRUN.
set -euo pipefail

: "${AGILE_DIR:?AGILE_DIR must be set}"
NP=${NP:-4}
ARCH=${ARCH:-gnu}
MPIRUN=${MPIRUN:-mpirun --oversubscribe}
OUTDIR=${OUTDIR:-../output/convergence}
SETUPS=${SETUPS:-"x diag"}

X_NLIST=${X_NLIST:-"50 100 200 400"}
X_LIMITERS=${X_LIMITERS:-"minmod vanleer mcbeta koren"}
X_STEPPERS=${X_STEPPERS:-"onestep twostep threestep"}
X_FLUXES=${X_FLUXES:-"tvdlf hll"}

D_NLIST=${D_NLIST:-"50 100 200"}
[ "${WITH400:-0}" = 1 ] && D_NLIST="$D_NLIST 400"
D_LIMITERS=${D_LIMITERS:-"vanleer koren"}
D_STEPPERS=${D_STEPPERS:-"threestep"}
D_FLUXES=${D_FLUXES:-"tvdlf"}

# AGILE needs an even block size >= 4. Per resolution, take the largest even
# divisor of N up to 50 that still leaves at least 2 blocks per direction:
# 50 -> 10, 100 -> 50, 200 -> 50, 400 -> 50 (large blocks keep the ghost-cell
# and geometry-array overhead down).
block_for () {
  local N=$1 b
  for b in 50 40 20 10 8 6 4; do
    if (( N % b == 0 && N / b >= 2 )); then echo $b; return; fi
  done
  echo "no even block size >= 4 divides N=$N into >= 2 blocks" >&2; exit 1
}

here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
mkdir -p "$OUTDIR"

echo "== building (arch=$ARCH)"
make arch="$ARCH" >/dev/null

run_one () {   # setup lim stepper flux N
  local setup=$1 lim=$2 stepper=$3 flux=$4 N=$5
  local name="conv_${setup}_${lim}_${stepper}_${flux}_N${N}"
  local base="$OUTDIR/$name"
  local cfl=0.7d0 np n2 n3 b2 b3 tend v nblocks maxb BLOCK
  BLOCK=$(block_for "$N")
  [ "$stepper" = onestep ] && cfl=0.4d0
  if [ "$setup" = x ]; then
    n2=10; n3=10; b2=10; b3=10; tend=1.0d0; v="1.0d0, 0.0d0, 0.0d0"; np=1
  else
    n2=$N; n3=$N; b2=$BLOCK; b3=$BLOCK; tend=0.333333333333333333d0
    v="1.0d0, 1.0d0, 1.0d0"; np=$NP
  fi
  nblocks=$(( (N/BLOCK) * (n2/b2) * (n3/b3) ))
  maxb=$(( (nblocks + np - 1)/np + 4 ))
  rm -f "$base.log"
  sed -e "s|@BASE@|$base|" -e "s|@TEND@|$tend|" -e "s|@STEPPER@|$stepper|" \
      -e "s|@FLUX@|$flux|" -e "s|@LIM@|$lim|" -e "s|@MAXB@|$maxb|" \
      -e "s|@B1@|$BLOCK|" -e "s|@B2@|$b2|" -e "s|@B3@|$b3|" \
      -e "s|@N1@|$N|" -e "s|@N2@|$n2|" -e "s|@N3@|$n3|" \
      -e "s|@CFL@|$cfl|" -e "s|@V@|$v|" -e "s|@SETUP@|$setup|" \
      conv_template.par > "$base.par"
  printf '%-48s ' "$name"
  local t0=$SECONDS
  if $MPIRUN -np "$np" ./agile -i "$base.par" > "$base.out" 2>&1; then
    printf 'done  (%ds)\n' $((SECONDS - t0))
  else
    printf 'FAILED, see %s.out\n' "$base"
  fi
  rm -f "$base"*.dat     # the end-state snapshot is not needed
}

for setup in $SETUPS; do
  if [ "$setup" = x ]; then
    lims=$X_LIMITERS; steps=$X_STEPPERS; fluxes=$X_FLUXES; nlist=$X_NLIST
  else
    lims=$D_LIMITERS; steps=$D_STEPPERS; fluxes=$D_FLUXES; nlist=$D_NLIST
  fi
  echo "== setup $setup"
  for lim in $lims; do for stepper in $steps; do for flux in $fluxes; do
    for N in $nlist; do run_one "$setup" "$lim" "$stepper" "$flux" "$N"; done
  done; done; done
done

python3 convergence.py "$OUTDIR"
