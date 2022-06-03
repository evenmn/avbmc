import avbmc as mc

system = mc.System()

# initialize from object
sampler = mc.Metropolis(system)
system.set_sampler(sampler)

# initialize by keyword (umbrella initializations should still work,
# albeit the function ntabulated and tabulated arguments are not used)
system.set_sampler("metropolis")
system.set_sampler("metropolis", lambda n : 0.1 * (n-20)**2)
system.set_sampler("metropolis", lambda n : 0.1 * (n-20)**2, 10)
system.set_sampler("metropolis", [0.1, 1.0])

