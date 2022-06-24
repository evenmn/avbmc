import avbmc as mc

# initialize box manually
system = mc.System(initialize=False)
box = mc.Box(system)
system.add_box(box)

# create ideal gas object
forcefield = mc.IdealGas(box, ["H"])
assert forcefield.ntype == 1
assert forcefield.label2type["H"] == 0
assert forcefield.unique_labels == ["H"]
box.set_forcefield(forcefield)
system.set_forcefield(forcefield)
system.set_forcefield(forcefield, box_id=0)

# initialize from keyword
system.set_forcefield("idealgas", ["H"])
