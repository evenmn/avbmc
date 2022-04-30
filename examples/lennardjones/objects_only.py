import avbmc as mc

system = mc.System()
box = mc.Box(system)
system.add_box(box)

forcefield = mc.LennardJones(box, "params.lj")
box.set_forcefield(forcefield)

box.add_particle("Ar", (1.,0.,0.))
box.snapshot("initial.xyz")

move = mc.Trans(system, box, 0.1)
system.add_move(move)
system.run_mc(100)
