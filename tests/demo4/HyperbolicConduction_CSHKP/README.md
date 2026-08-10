# Flare current-sheet transverse-conduction demonstration

This directory contains the user module, production namelists, and plotting
script for the 2.5D Yokoyama--Shibata flare current-sheet comparison used in
the AMRVAC 4.0 thermal-conduction section.  The modified user module is
`mod_usr.t`; it is paired with the fixed-resistivity configuration in
`production/flare_common.par` and `production/production_l3_t06.par`.

The four plotted cases are the matched parallel limit, electron-magnetization
closure, and effective ratios \(f_\perp=10^{-3}\) and \(10^{-2}\).  Their
overlays are `production/case00_parallel.par`,
`production/case01_electron_magnetization.par`,
`production/case03_effective_1e-3.par`, and
`production/case04_effective_1e-2.par`.

## Build, run, and plot

From this directory, set `AMRVAC_DIR` to the root of the AMRVAC checkout:

```bash
"$AMRVAC_DIR/setup.pl" -d=2 -v=3
make AMRVAC_DIR="$AMRVAC_DIR"
mkdir -p production/datamr
for case in \
  case00_parallel \
  case01_electron_magnetization \
  case03_effective_1e-3 \
  case04_effective_1e-2; do
  mpirun -np 4 ./amrvac -i production/flare_common.par \
    production/production_l3_t06.par "production/${case}.par"
done
python3 scripts/analyze_cshkp.py \
  --data-dir production/datamr --snapshot 5 --output-dir results
```

The analysis script reads the custom transport logs and the corresponding
AMRVAC snapshots, then writes CSV/NPZ source data and the five-panel figure to
`results/`.  The production calculation uses a \(128^2\) base mesh, AMR level
3, HLL fluxes, van Leer reconstruction, and the early-time comparison window
specified in the paper.  Python plotting additionally requires `numpy`,
`matplotlib`, and `yt`.
