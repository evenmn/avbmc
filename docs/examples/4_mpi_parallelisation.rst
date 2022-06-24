MPI parallelisation
===================

In a previous version of the code, MPI was used directly in the source files, making it possible to sample systems in parallel also when building from source without a Python interface. However, as this parallelisation method was found diffucult to combine with the Python interface, it was removed. Since Monte Carlo simulations are rather easy to parallelise, we instead encourage the users to parallelise the code directly in the Python script using :code:`mpi4py`. Below is a simple example where Lennard-Jones systems are sampled in parallel in the canonical ensemble. Run code on 4 processes using

.. code-block:: bash

   mpirun -n 4 python mpi_nvt.py


.. literalinclude:: ../../examples/lennardjones/mpi_nvt.py
   :language: python
