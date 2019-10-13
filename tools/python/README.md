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