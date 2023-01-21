
Adding constraints
------------------

The ``add_constraint`` function takes the following signature:

.. code-block:: python

   system.add_constraint(constraint, element1, element2, **kwargs)

where ``element1`` is the source element and ``element2`` is the target element. Example:

.. code-block:: python

   system.add_constraint("maxdistance", "O", "H", 2.0)

The example above ensures that there is always a particle of type ``H`` within a distance ``2.0`` from every particle of type ``O``. The logical counterpart would be

.. code-block:: python

   system.add_constraint("mindistance", "O", "H", 2.0)

The number of neighbors can also be restricted by the more general constraints ``MinNeigh`` and ``MaxNeigh``. To force every particle of type ``O`` to have two neighbors within a distance of 4.0 that are identical to itself, use

.. code-block:: python

   system.add_constraint("minneigh", "O", "O", 4.0, 2, box_id=0)

Dealing with constraint objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Like any other part of the code, the constraints are treated as objects that can be called directly in the Python script. For example, the ``MaxDistance`` object can be created by:

.. code-block:: python

   maxdistance = mc.MaxDistance(box, "O", "H", 2.0)

and then be set by:

.. code-block:: python

   box.add_constraint(maxdistance)

or

.. code-block:: python

   system.add_constraint(maxdistance, box_id=0)

The Stillinger constraint
^^^^^^^^^^^^^^^^^^^^^^^^^

Unlike most other constraints, the ``Stillinger`` constraint cannot be set using the soft ``add_constraint`` function presented above. The ``Stillinger`` constraint usually requires interactions with its methods, and should therefore always be manually created as an object. The pipeline would usually be

.. code-block:: python

   stillinger = mc.Stillinger(box)
   stillinger.set_criterion("O", "O", 4.0)
   stillinger.set_criterion("O", "H", 1.6)
   stillinger.set_criterion("H", "H", 0.0)
   box.add_constraint(stillinger)
