# junopy  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23898717.svg)](https://doi.org/10.5281/zenodo.23898717)
[![Documentation](https://img.shields.io/badge/docs-online-blue)](https://ysar.github.io/junopy/)

A python package to download, parse and plot particle and fields data from the NASA Planetary Data 
System for the Juno spacecraft.

# Installation

To install, clone this repository and install using `pip`.

```shell
git clone https://github.com/ysar/junopy.git --depth=1
cd junopy
```
And then,
```shell
pip install .
```
Or, using `uv` (recommended), 
```
uv pip install .
```
If you foresee making changes to `junopy`, using the `-e` (editable) flag - 
```
uv pip install -e .
```
# Usage

Please read the documentation in the `docs` directory or on the documentation [website](https://ysar.github.io/junopy). 
It will be particularly valuable to go through and run the examples. 

If you have `uv` installed, you can also build and serve the documentation website locally at 
[`http://127.0.0.1:8000`](http://127.0.0.1:8000) using (needs the `-e` flag during installation)

```shell
uv run mkdocs serve
```
This can also be done without using `uv`, but it involves installing `mkdocs`, `mkdocstrings`, and
`mkdocs-jupyter` into the environment.  

Using `junopy` is fairly straightforward with various imports (see example notebooks). One thing to 
remember is to use `junopy.config` to point towards the directory where Juno data is stored on your
disk. 

For example,
```python
from junopy import config
config.ConfigVars.set_root("../relative/path/to/some/folder/.junodata") 
```
This should be present in every script so that `junopy` can locate the previously downloaded SPICE 
kernels and other data. It is not necessary to have multiple `.junodata` folders on your 
system; it can be shared among different projects. User code, in the form of Jupyter notebooks 
or scripts, can be organized into different folders, all re-using the same `.junodata`.

> [!IMPORTANT]
> Verify that each script using `junopy` points to the correct `.junodata` directory.

If you are using a virtual environment (created either explicitly or implicitly using `uv`), you'll
need to create a "kernel" to link that environment to a Jupyter notebook. You only need to do this
once. Steps for doing this vary - see for example - [uv with Jupyter](https://docs.astral.sh/uv/guides/integration/jupyter/#creating-a-kernel).

# Capabilities
- **Magnetometer**: Any coordinate system, any cadence. Coordinate system conversions are seamless 
 if the relevant SPICE kernels are present.
- **Waves**: Survey-mode data can be analyzed (both E and B field). Burst-mode data may be added at
 a future date.
- **JADE**: Electron data processing is fairly mature. Ion SPECIES data can be processed in a similar
 manner to the electrons. TOF data can be read but processing is rudimentary.
- **JEDI**: Electron and ion TOFxE and TOFxPH data can be read and converted into energy-time and 
 pitch-angle time spectra. 
- Miscellaneous: Several other utilities are provided - calculating M-shell, making velocity-space 
 plots, etc.

# Contributing
Changes, bug fixes, and additions are very welcome. Please feel free to create an issue or pull request. 

# Disclaimer
This tool was built for personal research over several years and is not an official tool for NASA's 
Juno mission. All data used is publicly available from the PDS and NAIF. Errors and bugs may be 
present.

# Citation
Sarkango, Y. (2026). junopy: A Python package for processing particles and fields data from the Juno 
spacecraft [Computer software]. Zenodo. [https://doi.org/10.5281/zenodo.23898717](https://doi.org/10.5281/zenodo.23898717)
