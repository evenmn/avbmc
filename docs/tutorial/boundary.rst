
Setting boundary conditions
---------------------------

Every box is associated with boundary conditions, which can be set by the signature

.. code-block:: python

   system.set_boundary("boundary", box_id=-1, **kwargs)

By default, the boundaries are open, which may be set by

.. code-block:: python

   system.set_boundary("open", box_id=0)

Currently, ``Periodic`` is the only alternative to ``Open``\ , wraps the boundaries according to the minimum-image convention. The ``Periodic`` boundary requires box lengths to be defined:

.. code-block:: python

   system.set_boundary("periodic", length=[10, 10, 10], box_id=0)

The box will always have a corner in origin and spawn the space in positive directions. Make sure that the number of box dimensions match the number of dimensions of the particle positions.
