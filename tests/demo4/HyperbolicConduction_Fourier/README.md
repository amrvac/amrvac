# Oblique Fourier thermal-conduction demonstration

This directory is the reproducible input and analysis package for the Fourier
case used in the AMRVAC 4.0 thermal-conduction section.  It uses a periodic
two-dimensional unit square, a uniform field direction
\(\boldsymbol b=(2,1)/\sqrt{5}\), and equal-wavenumber parallel, mixed, and
perpendicular temperature modes.  The production matrix compares the
full-tensor RKL2 parabolic reference with the two-channel hyperbolic scheme at
three resolutions and three relaxation times.

The user module is `mod_usr.t`.  The base namelists are `common.par`,
`parabolic.par`, and `hyperbolic.par`; `scripts/generate_matrix.py` writes the
generated per-case overlays to `params/production/`.

## Build, run, and plot

From this directory, set `AMRVAC_DIR` to the root of the AMRVAC checkout:

```bash
"$AMRVAC_DIR/setup.pl" -d=2 -v=2
make AMRVAC_DIR="$AMRVAC_DIR"
python3 scripts/generate_matrix.py
python3 scripts/run_matrix.py --ranks 4
python3 scripts/analyze_benchmark.py
```

The last command reads `runs/*/output.log` and `output0000.dat`, writes source
tables under `results/source_data/`, and creates the SVG, PDF, PNG, and TIFF
figure under `results/figures/`.  The matrix runner stops on the first failed
case and does not overwrite an existing run directory; use `--resume` only
for completed cases.

The benchmark requires the AMRVAC build dependencies plus Python packages
`numpy`, `matplotlib`, and `yt`.  It is a uniform-grid verification case with
no AMR; the displayed values are dimensionless.
