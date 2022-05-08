import avbmc as mc
from tqdm import trange

system = mc.System()
system.set_temp(1.3)

system.set_forcefield("lennardjones", "params.lj")

system.add_particle("Ar", [0,0,0])

system.set_thermo(1, "mc.log", ["step", "poteng"])
system.set_dump(1, "mc.xyz", ["x", "y", "z"])

system.add_move("trans", dx=0.1)

# Customised iteration loop
system.initialize_mc_run()
for i in trange(100):
    system.run_mc_cycle();
