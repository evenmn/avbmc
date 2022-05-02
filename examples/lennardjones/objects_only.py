import avbmc as mc

system = mc.System(initialize=False)
box = mc.Box(system)
system.add_box(box)
system.set_temp(1.3)

forcefield = mc.LennardJones(box, "params.lj")
box.set_forcefield(forcefield)

box.add_particle("Ar", (1.,0.,0.))

box.set_thermo(1, "mc.log", ["step", "poteng"])
box.set_dump(1, "mc.xyz", ["x", "y", "z"])

move = mc.Trans(system, box, 0.1)
system.add_move(move)
system.run_mc(100)
