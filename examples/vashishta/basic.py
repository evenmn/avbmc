import avbmc as mc

system = mc.System()
system.set_temp(300.)

system.set_forcefield("vashishta", "H2O.wang.vashishta")
system.read_particles("water.xyz")

system.add_move("trans", dx=0.1)

system.set_dump(1000, "mc.xyz", ["x", "y", "z"])
system.set_thermo(1000, "mc.log", ["step", "atoms", "poteng"])

system.run_mc(10000000)
