
Adding particles
----------------

A single particle can be added to the box by its element type and position,

.. code-block:: python

   system.add_particle("H", [0, 0, 0])

For larger systems, it is recommended to initialize the system from an XYZ-file,

.. code-block:: python

   system.read_particles("particles.xyz")

or from a list of particle objects:

.. code-block:: python

   particles = mc.read_xyz("particles.xyz")
   system.add_particles(particles)

Make sure that all element types are covered by the forcefield parameter file - otherwise an error will be thrown. Again, the ``box_id``\ -argument might be used to specify which box to assign particles to, default is ``box_id=0``.

Initialize face-centered cube
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The positions of a face-centered cube can be written to a list by specifying the number of repeating unit cells in each direction, ``n``, and the length of a unit cell, ``d`` in the ``fcc`` function. Then, the particles can be added to the system by specifying the particle label. Currently, the face-centered cube initialization only support particles of the same type. An example:

.. code-block:: python

   positions = mc.fcc(n=3, d=1.7)
   system.add_particles("Ar", positions)
