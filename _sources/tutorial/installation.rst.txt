
Installing and importing
------------------------

Since this project is still experimentally, the library is not yet uploaded to PyPi. To install, run

.. code-block:: bash

   pip install git+https://github.com/evenmn/avbmc

or clone this repository and install from source by

.. code-block:: bash

   pip install .

This might take some time, as the entire library is built. When building is done, ``avbmc`` can be imported to Python:

.. code-block:: python

   import avbmc

or 

.. code-block:: python

   import avbmc as mc

In the examples later in the documentations, the latter line is assumed.

Initializing system
-------------------

To initialize a system, a new system object has to be created. It takes the following signature:

.. code-block:: python

   System(working_directory="", initialize=True)

The first argument, ``working_directory``\ , specifies the output file directory, while the second argument tells the code whether or not a box object should be initialized. Examples:

.. code-block:: python

   system = mc.System()  # initializing system with box writing to current directory
   system = mc.System("simulation")  # initializing system with box writing to "simulation" directory
   system = mc.System(initialize=False)  # initializing system with no box writing to current directory

A box can automatically be added to the system by

.. code-block:: python

   system.add_box()

or manually by

.. code-block:: python

   box = system.Box(system)
   system.add_box(box)

Every time a box is added, its ``box_id`` is increased by 1 compared to the previous initialized box, starting from 0. There is no upper limit for the number of boxes, but you are unlikely to need more than two.
