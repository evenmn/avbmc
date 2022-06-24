Grand-canonical ensemble studies with Vashishta
===============================================

The Vashishta potential is a class 3 forcefield, meaning that it allows for non-bonded three-body interactions. Thus, the potential is well-suited for the study of triatomic molecules. We will in this example show how water droplets can be sampled in the grand-canonical ensemble (:math:`\mu` PT) using AVBMC insertion and deletion moves of water molecules. We use umbrella sampling to sample areas of configuration space that one is usually unlikely to reach (in particular droplet sizes around maximum formation free energy). For practical reasons [...], the droplet is forced to be spherical by restricting all oxygen atoms to have at least three oxygen atoms as neighbors. 

.. literalinclude:: ../../examples/vashishta/upt.py
   :language: python
