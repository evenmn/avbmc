Grand-canonical simulations with Lennard-Jones
==============================================

Simulating a system in the grand-canonical ensemble (:math:`\mu` PT) is different from the canonical ensemble in the sense that we do not need fixed boundaries, but have to allow insertion and deletion moves. In this example, we use AVBMC insertion and deletion moves in addition to standard translational moves. We also restrict the system to one cluster to control the density.

.. literalinclude:: ../../examples/lennardjones/upt.py
   :language: python
