from setuptools import setup, find_packages
import os

project_name = "amrvac_pytools"
required = ['numpy', 'scipy', 'matplotlib', 'vtk']
version = '1.0.0'

setup(
    name = project_name,
    version = version,
    keywords = "interface data-analysis",
    python_requires = ">=3.6",
    install_requires = required,
    packages = find_packages(),
)