import avbmc as mc
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Set up system
system = mc.System()
system.set_temp(300.)
#system.set_chempot(-500.)
#system.set_seed(54456)

# Set forcefield and sampler
system.set_forcefield("vashishta", "params.vashishta")
system.set_sampler("umbrella", lambda n : 0.1 * (n-90)**2)

# Add constraints
system.add_constraint("minneigh", "O", "O", 3.0, 4)

# Read a single water molecule from file and initialize
water_mol = mc.read_xyz("water.xyz")
system.add_particles(water_mol)

# Add Monte Carlo moves
system.add_move("trans", 0.9, dx=0.1)
system.add_move("avbmcmol", 0.1, water_mol)

# Set outputs
system.set_dump(1000, f"dump_{rank}.xyz", ["x", "y", "z"])
#system.set_thermo(100, "mc.log", ["step", "atoms", "poteng"])

# Run Monte Carlo simulation
system.initialize_mc_run()
for i in range(200000):
    system.run_mc_cycle()

# Merge all histograms into one
size_histograms = comm.gather(system.get_size_histogram(), root=0)
if rank == 0:
    maxsize = 0
    for hist in size_histograms:
        if len(hist) > maxsize:
            maxsize = len(hist)

    cum_hist = np.zeros(maxsize, dtype=int)
    #for hist in size_histograms:
