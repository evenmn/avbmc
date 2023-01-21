System outputs
--------------

Writing box properties to file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Box properties can be written to a log file using the ``set_thermo`` method (inspired by LAMMPS' ``thermo`` command). The signature is the following:

.. code-block:: python

   system.set_thermo(freq, filename, outputs, box_id=0)

where ``outputs`` is a list of system properties, for instance ``step``\ , ``atoms``\ , ``poteng`` etc.. Example:

.. code-block:: python

   system.set_thermo(1, "mc1.log", ["step", "atoms", "poteng"], box_id=1)

The output file starts with a line of keywords, followed by a line of properties for every ``freq`` step. The file can easily be read using ```numpy.loadtxt`` <https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html>`_ or ```pandas.read_csv`` <https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html>`_.

Writing particle properties to file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similar to the system properties, particle properties can be written to a file with the signature

.. code-block:: python

   system.set_dump(freq, filename, outputs, box_id=0)

Possible particle properties is currently ``x``\ , ``y``\ , ``z``\ , ``xy`` and ``xyz``. An example:

.. code-block:: python

   system.set_dump(100, "mc0.xyz", ["xyz"], box_id=0)

The XYZ-file can be analysed using `Ovito <https://www.ovito.org/>`_\ , `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_\ , `ASE <https://wiki.fysik.dtu.dk/ase/>`_ or other software. NB: A low ``freq`` may slow down the run significantly. 

Writing number of particles histogram
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A histogram of system sizes throughout the run can be written to file with

.. code-block:: python

   system.write_size_histogram(filename, box_id=0)

To get the histogram as a Python list, use

.. code-block:: python

   histogram = system.get_size_histogram()

