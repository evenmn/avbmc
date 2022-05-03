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

