## Setting sampler
The signature of the `set_sampler` method is the following:
``` python
system.set_sampler(sampler, **kwargs)
```
The default sampler is `Metropolis`, which is set by
``` python
system.set_sampler("metropolis")
```
The only alternative that is currently available is the `Umbrella` sampler, which adds a weight to the system energy based on the number of particles in the system. This weight function is an arbitrary function $\mathbb{Z}$\Rightarrow$\mathbb{R}$, and is passed as an argument to the `Umbrella` constructor:
``` python
def f(n):
    return 0.003 * (n - 20)**2

system.set_sampler("umbrella", f=f, ntabulate=100)
```
where `ntabulate` is the maximum system size that is tabulated. The weight function can also be initialized by a `lambda`-function:
``` python
system.set_sampler("umbrella", f=lambda n : 0.003 * (n - 20)**2)
```

### Initialize umbrella from table (to be implemented)
In the future, we will also support umbrellas initialized from a table. A sketched, fictious example could look like:
``` python 
import numpy as np
x = np.linspace(-20, 100)
w = x**2 + 20

system.set_sampler("umbrella", w)
```
