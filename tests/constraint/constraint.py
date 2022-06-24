import avbmc as mc

# initialize system with ideal gas forcefield
system = mc.System(initialize=False)
box = mc.Box(system)
system.add_box(box)
system.set_forcefield("idealgas", ["H"])

# call box method
constraint = mc.MinNeigh(box, "H", "H", 3.0, 3)
box.add_constraint(constraint)
assert box.nconstraint == 1
assert len(box.constraints) == 1
box.rm_constraint(0)
assert box.nconstraint == 0
assert len(box.constraints) == 0

# call system method
system.add_constraint(constraint)
assert system.boxes[0].nconstraint == 1
assert len(system.boxes[0].constraints) == 1
system.rm_constraint(0)
assert system.boxes[0].nconstraint == 0
assert len(system.boxes[0].constraints) == 0

# call system method
system.add_constraint("minneigh", "H", "H", 3.0, 3)
assert system.boxes[0].nconstraint == 1
assert len(system.boxes[0].constraints) == 1
system.rm_constraint(0)
assert system.boxes[0].nconstraint == 0
assert len(system.boxes[0].constraints) == 0
