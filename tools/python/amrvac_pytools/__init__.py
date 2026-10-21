def _missing_dependency_loader(function_name, error):
    def _loader(*args, **kwargs):
        raise ImportError(
            "{} requires optional amrvac_pytools data-reading dependencies. "
            "Install tools/python with its full dependencies first."
            .format(function_name)
        ) from error

    return _loader


try:
    from .datfiles.reading.amrvac_reader import load_datfile
except ModuleNotFoundError as error:
    load_datfile = _missing_dependency_loader("load_datfile", error)

try:
    from .vtkfiles.read import load_vtkfile
except ModuleNotFoundError as error:
    load_vtkfile = _missing_dependency_loader("load_vtkfile", error)
