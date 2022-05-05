About AVBMC
========================
:code:`AVBMC` is a Python library for running Monte Carlo simulations of particles interacting with non-bonded potentials. 

Basic usage example for translational moves of two Argon atoms interacting with the Lennard-Jones potential.

.. code-block:: python 

    import avbmc as mc

    system = mc.System()

    system.set_temp(1.4)
    system.set_forcefield("lennardjones", "params.lj")

    system.add_particle("Ar", [0, 0, 0])
    system.add_particle("Ar", [1, 0, 0])
    system.add_move("trans", dx=0.1)

    system.run_mc(100)
