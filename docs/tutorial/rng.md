## Setting random number generator
The random number generator can be set by
``` python
system.set_rng("mt19937", seed=47383)
```
If the seed is not specified (or if it's negative), it will be given by a random device. Unfortunately, the Mersenne Twister pseudo random number generator is the only one implemented, and usually there is therefore no point of doing this. Contributions are welcome. To just set the seed, use
``` python
system.set_seed(47383)
```
or directly on the rng-object:
``` python
rng = mc.MersenneTwister;
rng.set_seed(47383)
system.set_rng(rng)
```
