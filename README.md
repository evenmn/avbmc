![build docs](https://github.com/henriasv/molecular-builder/workflows/build%20docs/badge.svg) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# AVBMC
Monte Carlo package for simulating systems of non-bonded particles. Open (uVT) and closed systems (NVT) are supported thought various move types, including translational moves and insertation/deletion moves (AVBMC). Additionally, a groups of particles can be inserted/removed, allowing non-bonded molecule to be inserted/removed. Code is written in C++.

## Prerequisites
The code depends on several features that were introduced in C++11. GCC 4.8.1 or any later version needed. MPI is needed for parallel processing. Nightly runs on macOS and Linux is performed using MPICH and OpenMPI.

## Build and run
There are two ways to build the code; with or without an input file parser. 

### With parser (recommended)
Build with:
```bash
make parser
```
and thereafter run with input script as an argument:
```bash
avbmc input.in
mpirun -n 4 avbmc input.in
```

### Without parser
Build with:
```bash
make dev
```
This procedure builds `src/main.cpp`, so this file has to be modified in order to configure the simulation. Thereafter, run with
```bash
make run
mpirun -n 4 ./main.out
```

## Configure simulations
The two build options above require different ways of configurating the simulations: With or without an input script.

### From input script
An input script might look like this:
```bash
# initialize box
set temp 0.7
set chempot -1.3
set boundary stillinger 1.5

# initialize positions
add particle Ar 0. 0. 0.
set mass Ar 1.
set forcefield lennardjones params.lj

# initialize Monte Carlo
set sampler umbrella square 32. 0.007
add move trans 0.94 0.1
add move avbmc 0.06 0.95 1.5

# set outputs
thermo 1 mc.log step atoms poteng acceptanceratio
dump 1 mc.xyz x y z

# run Monte Carlo simulation
take snapshot initial.xyz
run mc 100 1
take snapshot final.xyz
```

### Directly in main.cpp
The corresponding `main.cpp` looks like this:
```bash
#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include "box.h"
#include "forcefield/lennardjones.h"
#include "sampler/umbrella.h"
#include "moves/trans.h"
#include "moves/avbmc.h"

using namespace std;


int main()
{
    // initialize box
    Box box("simulation");
    box.set_temp(0.7);
    box.set_chempot(-1.3);

    // initialize particle at (0, 0, 0)
    box.add_particle("Ar", {0, 0, 0});
    box.set_mass("Ar", 1.);
    box.set_forcefield(new LennardJones(&box, "params.lj"));

    // initialize Monte Carlo
    double k = 0.007;
    double nc = 32.;
    auto f = [nc, k] (const int n) { return (k * (n - nc) * (n - nc)); };
    box.set_sampler(new Umbrella(&box, f));
    box.add_move(new Trans(&box, 0.01), 0.94);
    box.add_move(new AVBMC(&box, 0.9, 1.5), 0.06);

    // set outputs
    box.set_dump(1, "mc.xyz", {"x", "y", "z"});
    box.set_thermo(1, "mc.log", {"step", "atoms", "poteng", "acceptanceratio"});

    // run Monte Carlo simulation
    box.snapshot("initial.xyz");
    box.run_mc(10000, 1);
    box.snapshot("final.xyz");

    return 0;
}
```

## File formats
We try to stick to standard input file formats. For datafiles, the xyz-format is supported both to initialize the system and for dump files. The parameter files follow the LAMMPS standard. Both xyz-files and LAMMPS parameter files are oriented around chemical symbols (in contrast to for instance LAMMPS data files, which are oriented around particle types). This means that the user never needs to deal with the particle type indices, only the chemical symbols. Particles with the same chemical symbol will have the same properties (mass, charge and other properties given by the parameter file).

## License
