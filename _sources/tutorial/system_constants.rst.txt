
Setting system constants
------------------------

The code is currently able to simulate systems in the NVT and uVT ensembles. Constant temperature is set by 

.. code-block:: python

   system.set_temp(300.)

If uVT ensemble is used, the chemical potential is set by

.. code-block:: python

   system.set_chempot(10.)

Please be aware of the units - The units used here need to match the units used in the potential parameter file.
