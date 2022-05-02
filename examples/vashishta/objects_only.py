import avbmc as mc

system = mc.System(initialize=False)
box = mc.Box(system)
system.add_box(box)
system.set_temp(300.)

forcefield = mc.Vashishta(box, "H2O.wang.vashishta")
box.set_forcefield(forcefield)

box.read_particles("water.xyz")

move = mc.Trans(system, box, 0.1)
system.add_move(move)

box.set_dump(1000, "mc.xyz", ["x", "y", "z"])
box.set_thermo(1000, "mc.log", ["step", "atoms", "poteng"])

system.run_mc(10000000)
