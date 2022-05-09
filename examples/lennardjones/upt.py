import avbmc as mc

system = mc.System()
system.set_temp(1.3)
system.set_chempot(10.)

system.set_forcefield("lennardjones", "params.lj")

system.add_particle("Ar", [1, 0, 0])

system.add_move("trans", prob=0.5, dx=0.1)
system.add_move("avbmc", prob=0.5, particle="Ar", r_below=0.95, r_above=3.0)

system.run_mc(100)
