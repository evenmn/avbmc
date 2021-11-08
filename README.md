# AVBMC
Monte Carlo package for estimating the free energy and surface tension of water droplets of various sizes. Code is written in C++.

## Input files
We try to stick to standard input file formats. For datafiles, the xyz-format is supported both to initialize the system and for dump files. The parameter files follow the LAMMPS standard. Both xyz-files and LAMMPS parameter files are oriented around chemical symbols (in contrast to for instance LAMMPS data files, which are oriented around particle types). This means that the user never needs to deal with the particle type indices, only the chemical symbols. Particles with the same chemical symbol will have the same properties (mass, charge and other properties given by the parameter file).
