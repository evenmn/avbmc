[![build](https://github.com/evenmn/avbmc/actions/workflows/c-cpp.yml/badge.svg)] [![build docs](https://github.com/evenmn/avbmc/actions/workflows/docs.yml/badge.svg)] [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# AVBMC
Python library for non-bonded Monte Carlo simulations. Open (uVT) and closed systems (NVT) are supported thought various move types, including translational moves and insertation/deletion moves (AVBMC). Additionally, a groups of particles can be inserted/removed, allowing non-bonded molecule to be inserted/removed. The library is entirely written in C++ for performance reasons, and binded to Python using [pybind11](https://pybind11.readthedocs.io/en/stable/index.html).

## Prerequisites
The code depends on several features that were introduced in C++11. GCC 4.8.1 or any later version needed. Python 3.5 or later is required.

## Installation
As this library is still under development and can be seen as experimental, it is not yet uploaded to PyPi. Meanwhile, install from source using
``` bash
pip install git+https://github.com/evenmn/avbmc
```
It might takes some time, as the entire library is built.

## Basic usage
``` python
import avbmc as mc

system = mc.System()

system.set_temp(1.4)
system.set_forcefield("lennardjones", "params.lj")

system.add_particle("Ar", [0, 0, 0])
system.add_particle("Ar", [1, 0, 0])
system.add_move("trans", dx=0.1)

system.run_mc(100)
```
For more example usage, see the [documentation](https://evenmn.github.io/avbmc).

## File formats
We try to stick to standard input file formats. For datafiles, the xyz-format is supported both to initialize the system and for dump files. The parameter files follow the LAMMPS standard. Both xyz-files and LAMMPS parameter files are oriented around chemical symbols (in contrast to for instance LAMMPS data files, which are oriented around particle types). This means that the user never needs to deal with the particle type indices, only the chemical symbols. Particles with the same chemical symbol will have the same properties (mass, charge and other properties given by the parameter file).

## License
