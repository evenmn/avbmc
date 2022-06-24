import avbmc as mc

# create system
system = mc.System(initialize=False)

# add box to system
system.add_box()
assert system.nbox == 1

# create box object and add to system
box = mc.Box(system)
assert system.nbox == 1
assert len(system.boxes) == 1
system.add_box(box)
assert system.nbox == 2
assert len(system.boxes) == 2

# remove boxes
system.rm_box(1)
assert system.nbox == 1
assert len(system.boxes) == 1
system.rm_box(0)
assert system.nbox == 0
assert len(system.boxes) == 0
