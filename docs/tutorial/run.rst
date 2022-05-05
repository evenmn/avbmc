
Running the code
-------------------------------

Monte Carlo simulations can be run by calling the function :code:`run_mc`:

.. code-block:: python

   system.run_mc(cycles, move_per_cycle=1)

Making your own sampling loop (to be implemented)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   sampler = mc.Metropolis()
   system.set_sampler(sampler)

   for step in range(100):
       sampler.run(1)

.. code-block:: python

   from tqdm import tqdm

   for step in tqdm(range(100)):
       sampler.run(1)


Running molecular dynamics simulations (to be implemented)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   system.run_md(steps)
