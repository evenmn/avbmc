import avbmc as mc

system = mc.System()
system.set_temp(4.)

system.set_forcefield("lennardjones", "params.lj")
system.set_boundary("periodic", [6.8, 6.8, 6.8])

positions = mc.fcc(n=4, L=6.8)
system.add_particles("Ar", positions)
system.snapshot("initial.xyz")

system.set_dump(100, "mc.xyz", ["x", "y", "z"])

system.add_move("trans", dx=0.1)
system.run_mc(10000)
