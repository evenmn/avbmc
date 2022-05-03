# General documentation

## Installing and importing
Since this project is still experimentally, the library is not yet uploaded to PyPi. To install, run
``` bash
pip install git+https://github.com/evenmn/avbmc
```
or clone this repository and install from source by
``` bash
pip install .
```
This might take some time, as the entire library is built. When building is done, `avbmc` can be imported to Python:
``` python
import avbmc
```
or 
``` python
import avbmc as mc
```
In the examples later in the documentations, the latter line is assumed.

## Initializing system
To initialize a system, a new system object has to be created. It takes the following signature:
``` python
System(working_directory="", initialize=True)
```
The first argument, `working_directory`, specifies the output file directory, while the second argument tells the code whether or not a box object should be initialized. Examples:
``` python
system = mc.System()  # initializing system with box writing to current directory
system = mc.System("simulation")  # initializing system with box writing to "simulation" directory
system = mc.System(initialize=False)  # initializing system with no box writing to current directory
```
A box can automatically be added to the system by
``` python
system.add_box()
```
or manually by
``` python
box = system.Box(system)
system.add_box(box)
```
Every time a box is added, its `box_id` is increased by 1 compared to the previous initialized box, starting from 0. There is no upper limit for the number of boxes, but you are unlikely to need more than two.

## Setting system constants
The code is currently able to simulate systems in the NVT and uVT ensembles. Constant temperature is set by 
``` python
system.set_temp(300.)
```
If uVT ensemble is used, the chemical potential is set by
``` python
system.set_chempot(10.)
```
Please be aware of the units - The units used here need to match the units used in the potential parameter file.

## Setting forcefield
The forcefield is usually set by the signature
``` python
system.set_forcefield(forcefield_name, parameter_file)
```
If several boxes are found, the same forcefield is set to all boxes. Avoid this by using the `box_id`-argument. Example:
``` python
system.set_forcefield("lennardjones", "params.lj", box_id=1)
```
The `box_id` starts from 0 for the first box and increases with 1 for every box added.

### Dealing with forcefield-object
Since everything is threated as objects under the hood, it is also possible to deal with the raw objects:
``` python
forcefield = mc.LennardJones(box, "params.lj")
box.set_forcefield(forcefield)
```
### Ideal gas
An exception from the above examples is the ideal gas forcefield, which does not take a parameter file. It is still necessary to specify all the particle types in the box. Set ideal gas by
``` python
system.set_forcefield("idealgas", ["Si", "O"])
```

## Adding particles
A single particle can be added to the system by its element type and position,
``` python
system.add_particle("H", [0, 0, 0])
```
For larger systems, it is recommended to initialize the system from an XYZ-file,
``` python
system.read_particles("particles.xyz")
```
or from a list of particle objects:
``` python
particles = mc.read_xyz("particles.xyz")
system.add_particles(particles)
```
Make sure that all element types are covered by the forcefield parameter file - otherwise an error will be thrown. Again, the `box_id`-argument might be used to specify which box to assign particles to, default is `box_id=0`.

### Initialize face-centered cube
TO BE ADDED

## Adding moves
The signature of the `add_move` method is as follows:
``` python
system.add_move(move_name, prob=1.0, box_id=-1, **kwargs)
```
where `box_id=-1` means that the move is assiged to all available boxes. Some examples are:
``` python
system.add_move("trans", prob=0.5, dx=0.1)
system.add_move("avbmc", prob=0.25, particle="Si", r_below=0.9, r_above=3.0, box_id=0)
system.add_move("avbmcmol", prob=0.25, molecule=water_mol, r_below=0.95, r_above=3.0, r_inner=1.3)
```

## System outputs

### Writing system properties to file
System properties can be written to a log file using the `set_thermo` method (inspired by LAMMPS' `thermo` command). The signature is the following:
``` python
system.set_thermo(freq, filename, outputs, box_id=0)
```
where `outputs` is a list of system properties, for instance `step`, `natom`, `poteng` etc.. Example:
``` python
system.set_thermo(1, "mc1.log", ["step", "natom", "poteng"], box_id=1)
```
The output file starts with a line of keywords, followed by a line of properties for every `freq` step. The file can easily be read using [`numpy.loadtxt`](https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html) or [`pandas.read_csv`](https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html).

### Writing particle properties to file
Similar to the system properties, particle properties can be written to a file with the signature
``` python
system.set_dump(freq, filename, outputs, box_id=0)
```
Possible particle properties is currently `x`, `y`, `z`, `xy` and `xyz`. An example:
``` python
system.set_dump(100, "mc0.xyz", ["xyz"], box_id=0)
```
The XYZ-file can be analysed using [Ovito](https://www.ovito.org/), [VMD](https://www.ks.uiuc.edu/Research/vmd/), [ASE](https://wiki.fysik.dtu.dk/ase/) or other software. NB: A low `freq` may slow down the run significantly. 

### Writing number of particles histogram
TO BE WRITTEN

## Adding constraints
TO BE WRITTEN

## Setting sampler
TO BE WRITTEN

## Setting random number generator
TO BE WRITTEN
