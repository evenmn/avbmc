import avbmc as mc
import matplotlib.pyplot as plt

# Set up system
system = mc.System()
system.set_temp(300.)
#system.set_chempot(-50.)
#system.set_seed(54456)

# Set forcefield and sampler
system.set_forcefield("vashishta", "params.vashishta")
system.set_sampler("umbrella", lambda n : 0.1 * (n-20)**2)

# Read a single water molecule from file and initialize
water_mol = mc.read_xyz("water.xyz")
system.add_particles(water_mol)

# Add Monte Carlo moves
system.add_move("trans", 0.9, dx=0.1)
system.add_move("avbmcmol", 0.1, water_mol)

# Add constraints
#system.add_constraint("minneigh", "O", "O", 3.0, 4)

# Set outputs
system.set_dump(1000, "dump.xyz", ["x", "y", "z"])
#system.set_thermo(100, "mc.log", ["step", "atoms", "poteng"])

system.run_mc(500000)
"""
# summary
print("#ndrawn   #naccept   time")
for move in system.moves:
    print(move, move.ndrawn, move.naccept, move.cum_time)

# plot size histogram
system.write_size_histogram("size_histogram.txt")
size_histogram = system.get_size_histogram()

plt.plot(size_histogram[::3])
plt.show()
"""
