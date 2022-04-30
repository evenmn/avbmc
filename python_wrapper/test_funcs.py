import avbmc_core as avbmc

# initialize system
system = avbmc.System("dir")
system.set_temp(320.)
system.set_chempot(10.)
rng = avbmc.MersenneTwister()
system.set_rng(rng)
sampler = avbmc.Metropolis(system)
system.set_sampler(sampler)

# initialize box
box = avbmc.Box(system, 2)
boundary = avbmc.Open(box)
box.set_boundary(boundary)
forcefield = avbmc.Vashishta(box, "H2O.nordhagen.vashishta")
box.set_forcefield(forcefield)
box.read_particles("water.xyz")
system.add_box(box)

# set constraints
stillinger = avbmc.Stillinger(box, 1.0)
stillinger.set_criterion("O", "O", 3.0)
stillinger.set_criterion("O", "H", 1.6)
stillinger.set_criterion("H", "H", 0.0)
box.add_constraint(stillinger)
minneigh = avbmc.MinNeigh(box, "O", "O", 3.0, 2)
box.add_constraint(minneigh)

# add moves
move = avbmc.Trans(system, box, 0.1)
system.add_move(move, 1.0)

# set output
#box.set_dump(1, "mc.xyz", ["x", "y", "z"], True)
#box.set_thermo(1, "mc.log", ["step", "atoms", "poteng"], True)

# run Monte Carlo simulation
box.snapshot("initial.xyz", True)
print(system.temp)
system.run_mc(10000, 1)
print("her")
