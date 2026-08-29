# amrvac_pytools
## Installation

This python package is regularly updated, so we recommend installing it with "develop" arguments so that it does not need to be reinstalled with each update. This can be achieved in 3 distinct ways
> *nb:* all methods require `cd amrvac_pytools-project` is executed beforehand

- **recommended:** using `conda`
```bash
conda install conda-build
conda develop .
```
- using `pip`
```bash
pip install -e .
```
- or simply vanilla python:
```bash
python setup.py develop
```

## usage

```python
import amrvac_pytools
```

### load data
```python
df = amrvac_pytools.load_vtkfile("path/to/vtk/file")
# or
df = amrvac_pytools.load_datfile("path/to/dat/file")
```

### writeme

## Data-driven boundary preprocessing

The `amrvac_pytools.datadriven` package provides Python preprocessing and case
staging for AMRVAC data-driven simulations. Two notebook entry points separate
the single-frame and time-dependent workflows:

```bash
jupyter lab tools/python/notebooks/DataConstrain.ipynb
jupyter lab tools/python/notebooks/DataDriven.ipynb
```

Both accept raw HMI `field/inclination/azimuth/disambig` data or prepared
SHARP/CEA `Br/Bt/Bp` data. `DataConstrain.ipynb` stages the single-frame chain:

```text
selected frame -> optional shared preprocessing -> optional Grad-Rubin alpha cleaning
-> PotentialField -> one selected NLFFF method -> DataConstrained
```

`DataDriven.ipynb` operates on a time sequence. It uses the first selected
frame as the reference field, stages the same potential-plus-one-NLFFF initial
field, and then stages one of two B-only evolution modes:

```text
magnetic sequence -> optional shared preprocessing -> optional Grad-Rubin alpha cleaning
-> PotentialField -> one selected NLFFF method -> TMF or DataDriven MHD
```

The public notebooks expose exactly one post-potential NLFFF selection:
`legacy_mfr`, `optimization`, or `grad_rubin`. Shared vector preprocessing is
controlled only by `PREPROCESS_VECTOR=False` by default; it is optional for the
legacy embedded-MHD MFR path and recommended for observational Optimization and
Grad-Rubin runs. Grad-Rubin alpha cleaning is controlled only by
`GR_ALPHA_CLEANING=False` by default; when enabled, the external-alpha product
is generated from the same already-preprocessed boundary field and is marked
unused if another NLFFF method is selected.

All three NLFFF methods expose the same common `<base_filename>_nlfff_metrics.csv`
path by default. Method-specific detailed CSV histories such as legacy
`_mflog.csv` remain opt-in through explicit detailed-history switches; older
legacy MFR logs can be normalized into the common metrics schema with unavailable
quantities recorded as `NaN` and documented in an adapter audit JSON.

The V1 time-dependent workflow writes one `B_XXXX.dat` file per observation.
It does not write velocity files and does not require DAVE/DAVE4VM. Observation
times remain in the boundary files; any time acceleration is an AMRVAC case
parameter.

The notebook uses a small orchestration API. Lower-level preparation and case
staging functions remain public for scripts and future time-dependent workflows.

```python
from amrvac_pytools.datadriven import (
    create_data_constrain_workflow,
    create_data_driven_workflow,
    stage_grid_config,
)

workflow = create_data_constrain_workflow(
    amrvac_root=AMRVAC_ROOT,
    project_dir=PROJECT_DIR,
    input_dir=INPUT_DIR,
    relaxation_grid=stage_grid_config(2, 1),
    evolution_grid=stage_grid_config(1, 2),
)

time_dependent = create_data_driven_workflow(
    amrvac_root=AMRVAC_ROOT,
    project_dir=PROJECT_DIR,
    input_dir=INPUT_DIR,
    relaxation_grid=stage_grid_config(2, 1),
    evolution_grid=stage_grid_config(1, 2),
)
```

The package is split by responsibility: `pipeline.py` prepares magnetic data,
`cases.py` stages AMRVAC cases and parameter files, `diagnostics.py` analyzes
magnetofrictional convergence, and `workflow.py` plus
`data_driven_workflow.py` provide notebook-facing orchestration. Sequence APIs
write one `B_XXXX.dat` file per frame. Each project also contains `project.json`,
which records shared input and grid configuration plus independent
DataConstrain and DataDriven workflow state for reuse after moving the project.

Velocity estimation is intentionally an extension point in V1. DAVE/DAVE4VM is
not vendored and is not imported by default; B-V driving can be added later via
the `VelocityEstimator` adapter interface.
