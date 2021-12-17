import numpy as np
import matplotlib.pyplot as plt

filename = "mc.log"

data = np.loadtxt(filename)

# Plot number of particles
plt.plot(data[:, 0], data[:, 1])
plt.show()

# Number histogram
max_ = int(np.max(data[:, 1]))
n, bins, patches = plt.hist(data[:, 1], max_)
plt.show()

# dG histogram
def w(n):
    return 0.012 * (n - 32) ** 2

n_ = np.asarray(range(1, max_+1))
dG = 0.7 * np.log(n[0]/n)
plt.plot(n_, dG, label="apparent dG")
plt.plot(n_, w(n_), label="umbrella")
plt.plot(n_, dG - w(n_), label="dG")
plt.legend(loc='best')
plt.show()
