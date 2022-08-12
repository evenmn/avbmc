Statistics
----------

Sampling statistics can tell a lot about the simulation and whether or not it was successful. One of the most important metrics for this is the acceptance ratio, which can indicate if the step length was appropriate. In the code, the number of attempts and accepted attempts are stored for each move. They can be found by :code:`move.ndrawn` and :code:`move.naccept`, respectively. Additionally, the cummulative CPU-time for each move is stored in :code:`move.cum_time`. The move statistics can be output by

.. code-block:: python

   print("Label #drawn #accepted acceptance time")
   for move in system.moves:
       print(move, move.ndrawn, move.naccept, move.naccept/move.ndrawn, move.cum_time)

Not much more statistics is stored in the code itself, but data can easily be stored in Python by using a customized sampling loop. For instance, if one wants to use resampling to estimate the variance of the energy, one would need to store energies after all cycles. This can be done by

.. code-block:: python

   potengs = []
   for cycle in range(1000):
       system.run_mc_cycle()
       potengs.append(system.poteng)

Other variables can be stored similarly throughout the sampling.

