import avbmc as mc

system = mc.System()

# initialize from function
def f(n):
    return 0.1 * (n-20)**2
sampler = mc.Umbrella(system, f, ntabulated=100)
system.set_sampler(sampler)

# initialize from lambda function
sampler = mc.Umbrella(system, lambda n : 0.1 * (n-20)**2, ntabulated=100)
system.set_sampler(sampler)

# initialize from tabulated elements
tabulated = [f(n) for n in range(100)]
sampler = mc.Umbrella(system, tabulated)
system.set_sampler(sampler)

# initialize by keyword
system.set_sampler("umbrella", f)
system.set_sampler("umbrella", f, ntabulated=100)
system.set_sampler("umbrella", tabulated)

