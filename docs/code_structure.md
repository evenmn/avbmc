# Code structure
The AVBMC code contains several thousand lines of code, and no matter 
how much effort we put in making it clean, it still can easily be 
confusing. It is oriented about the `System` class, which manages 
variables that are global for the entire system. This includes the
dimensionality, different moves, the sampler and random number generator.
Also, the system consists of one or several boxes.

Each box is initialized with a boundary object giving the boundary 
conditions of the box. The current boundary conditions are `Periodic`
and `Open`. Furthermore, all particle types in the box have to be covered by the
`ForceField` object. `LennardJones`, `Vashishta` and `IdealGas` are the
available force-fields. Particle constraints can be assigned to the box from a
`Constraint` object. Implemented constraints include `MaxDistance`,
`MinDistance`, `MaxNeigh`, `MinNeigh` and `Stillinger`. For molecular dynamics
runs, the box integrator has initialized from the `Integrator` set of classes, 
including `Euler`, `EulerCromer`, `VelocityVerlet` and `RungeKutta4`. 

For Monte Carlo runs, one or more moves have to be specified from the `Moves`
class. The `Moves` class belongs to the system rather than a box to make
inter-box moves possible. The available moves are `Trans`, `TransMH`, `AVBMCIn`,
`AVBMCOut`, `AVBMCMolIn` and `AVBMCMolOut`. 

The sampling technique also have to be defined before a Monte Carlo run. Code 
supports Metropolis sampling and umbrella sampling through the classes
`Umbrella` and `Metropolis`. 

## Code Chart
```
└── System
    ├── Box
    │   ├── Particle
    │   ├── DistanceManager
    │   ├── Boundary
    │   │   ├── Open
    │   │   ├── Periodic
    │   │   ...
    │   ├── ForceField
    │   │   ├── LennardJones
    │   │   ├── Vashishta
    │   │   ...
    │   ├── Integrator
    │   │   ├── Euler
    │   │   ├── VelocityVerlet
    │   │   ...
    │   └── Constraint
    │       ├── Stillinger
    │       ├── MinNeigh
    │       ...
    ├── Moves
    │   ├── Trans
    │   ├── AVBMCIn
    │   ...
    ├── Sampler
    │   ├── Metropolis
    │   ├── Umbrella
    │   ...
    └── RandomNumberGenerator
        ├── MersenneTwister
        ...
```

## Extending Library
The code has be written to make it easy to extend it without changing the basic
framework. This including adding new boundary conditions, force-fields,
integrators, constraints, moves, samplers and random number generators. To add
a class, look at the existing classes and try to modify it for your purposes.

## Changing Framework
If the current framework does not match your needs, please contact the
developers.

## Other Classes
There is a bunch of classes not mentioned yet, that is `Thermo` and `Dump`, and
a couples of utility files `init_position.cpp` and `init_velocity.cpp`.
