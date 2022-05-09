import avbmc as mc
from mpi4py import MPI
import matplotlib.pyplot as plt

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Set up system
system = mc.System()
system.set_temp(300.)

system.set_forcefield("lennardjones", "params.lj")
system.set_boundary("periodic", [6.8, 6.8, 6.8])

positions = mc.fcc(n=4, L=6.8)
system.add_particles("Ar", positions)
if rank == 0:
    system.snapshot("initial_fcc.xyz")

system.set_dump(100, f"dump_{rank}.xyz", ["x", "y", "z"])
system.set_thermo(1, f"log_{rank}.txt", ["step", "poteng"])

system.add_move("trans", dx=0.1)
system.initialize_mc_run()
for i in range(100000):
    system.run_mc_cycle()
