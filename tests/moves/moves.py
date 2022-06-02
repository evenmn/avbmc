import avbmc as mc

system = mc.System()

system.set_forcefield("idealgas", ["H"])

system.add_move("trans")
assert system.nmove == 1
system.rm_move(0)
assert system.nmove == 0
