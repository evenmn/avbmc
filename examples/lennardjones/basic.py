import avbmc as mc

system = mc.System()
system.set_temp(1.3)

system.set_forcefield("lennardjones", "params.lj")

system.add_particle("Ar", [0,0,0])

system.set_thermo(1, "mc.log", ["step", "poteng"])
system.set_dump(1, "mc.xyz", ["x", "y", "z"])

system.add_move("trans", dx=0.1)
system.run_mc(100)
