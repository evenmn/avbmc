Inter-box swap moves
====================

(To be implemented)

import avbmc as mc

system = mc.System()
system.set_seed(45820)
system.set_temp(2.5)

system.set_forcefield("lennardjones", "params.lj")
positions = mc.fcc(4, 6.8)
system.add_particles("Ar", positions)
system.set_boundary("periodic", [6.8, 6.8, 6.8])

# thermalize
system.add_move("trans", dx=0.5)
system.run_mc(100000)
system.rm_moves()
#system.rm_move(0)

# add another box
system.add_box()
system.set_forcefield("idealgas", ["Ar"], box_id=1)
system.add_particles("Ar", positions, box_id=1)

system.add_move("avbmcswap", 0.5, "Ar", box_id=0, box_id2=1)
system.add_move("trans", 0.5, dx=0.5, box_id=0)
system.run_mc(100000)
