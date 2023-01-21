
Setting forcefield
------------------

The forcefield is usually set by the signature

.. code-block:: python

   system.set_forcefield(forcefield_name, parameter_file)

If several boxes are found, the same forcefield is set to all boxes. Avoid this by using the ``box_id``\ -argument. Example:

.. code-block:: python

   system.set_forcefield("lennardjones", "params.lj", box_id=1)

The ``box_id`` starts from 0 for the first box and increases with 1 for every box added.

Dealing with forcefield-object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since everything is treated as objects under the hood, it is also possible to deal with the raw objects:

.. code-block:: python

   forcefield = mc.LennardJones(box, "params.lj")
   box.set_forcefield(forcefield)

Ideal gas
^^^^^^^^^

An exception from the above examples is the ideal gas forcefield, which does not take a parameter file. It is still necessary to specify all the particle types in the box. Set ideal gas by

.. code-block:: python

   system.set_forcefield("idealgas", ["Si", "O"])
