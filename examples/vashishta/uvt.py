import avbmc as mc

system = mc.System()
system.set_temp(300.)

system.set_forcefield("vashishta", "H2O.wang.vashishta")

water_mol = mc.read_xyz("water.xyz")
system.add_particles(water_mol)

system.add_move("trans", 0.9, dx=0.1)
system.add_move("avbmcmol", 0.5, water_mol)

system.set_dump(100, "mc.xyz", ["x", "y", "z"])
system.set_thermo(100, "mc.log", ["step", "atoms", "poteng"])

system.run_mc(10000)
