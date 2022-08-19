Statistics
----------

Sampling statistics can tell a lot about the simulation and whether or not it was successful. One of the most important metrics for this is the acceptance ratio, which can indicate if the step length was appropriate. In the code, the number of attempts and accepted attempts are stored for each move. They can be found by :code:`move.ndrawn` and :code:`move.naccept`, respectively. Additionally, the cummulative CPU-time for each move is stored in :code:`move.cum_time`. The move statistics can be output by

.. code-block:: python

   print("Label #drawn #accepted acceptance time")
   for move in system.moves:
       print(move, move.ndrawn, move.naccept, move.naccept/move.ndrawn, move.cum_time)

On can also call the built-in function :code:`system.print_statistics()`, which outputs a neat table to the terminal. It might look like this:

.. code-block:: bash

   ------------------------------------------------------------------------
   |    Move     | #drawn | #accept | #reject | acc. ratio | CPU-time (s) |
   ------------------------------------------------------------------------
   | Trans       |    492 |     120 |     372 |   0.243902 |     0.711072 |
   | AVBMCMolIn  |    260 |     104 |     156 |   0.400000 |     0.076163 |
   | AVBMCMolOut |    249 |       0 |     249 |   0.000000 |     0.119231 |
   ------------------------------------------------------------------------

For a deeper analysis, more columns can be specified with an argument list:

.. code-block:: python

   system.print_statistics(["move", "ndrawn", "nreject", "rejtarg", "rejout"])

which might output something like this:

.. code-block:: bash

   -----------------------------------------------------------------
   |    Move     | #drawn | #reject | #reject out | #reject target |
   -----------------------------------------------------------------
   | Trans       |    498 |     382 |           - |              - |
   | AVBMCMolIn  |    235 |     150 |           - |              0 |
   | AVBMCMolOut |    268 |     268 |           0 |              0 |
   -----------------------------------------------------------------

The possible columns are :code:`move`, :code:`ndrawn`, :code:`naccept`, :code:`nreject`, :code:`accratio`, :code:`cputime`, :code:`rejout` and :code:`rejtarg`.

A LaTeX style is also under development, which is expected to work like this:

.. code-block:: python

   latex_table = system.print_statistics(style="latex", print=false)


Custom statistics
^^^^^^^^^^^^^^^^^
Not much more statistics is stored in the code itself, but data can easily be stored in Python by using a customized sampling loop. For instance, if one wants to use resampling to estimate the variance of the energy, one would need to store energies after all cycles. This can be done by

.. code-block:: python

   potengs = []
   for cycle in range(1000):
       system.run_mc_cycle()
       potengs.append(system.poteng)

Other variables can be stored similarly throughout the sampling.

