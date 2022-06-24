import avbmc as mc

system = mc.System()

# create rng object
rng = mc.MersenneTwister(seed=12345)
rng = mc.MersenneTwister(seed=-1)
rng = mc.MersenneTwister()

# test RNG
N = 1000
probabilities1 = [1.]
probabilities2 = [0., 0., 1.]
probabilities3 = [0., 0.5, 0., 0.5, 0.]
to_shuffle = [1, 2, 3, 4, 5]
for _ in range(N):
    i = rng.next_int(10)
    assert 0 <= i < 10
    d = rng.next_double()
    assert 0 < d < 1
    d = rng.next_gaussian(0, 0)
    assert (d - 0.0) < 1e-14
    i = rng.choice(probabilities1)
    assert i == 0
    i = rng.choice(probabilities2)
    assert i == 2
    i = rng.choice(probabilities3)
    assert i == 1 or i == 3
    shuffled = rng.shuffle(to_shuffle)
    assert len(shuffled) == len(to_shuffle)

# initialize from object
system.set_rng(rng)
    
# initialize by keyword
system.set_rng("mersennetwister")
system.set_rng("mt19937")

# set seed
system.set_seed(12345)
