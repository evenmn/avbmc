import avbmc as mc

# initialize box manually
system = mc.System(initialize=False)
box = mc.Box(system)
system.add_box(box)
system.set_forcefield("idealgas", ["H"])

# add and remove particle object
particle = mc.Particle("H", [0,0,0])
box.add_particle(particle)
assert box.npar == 1
box.rm_particle(0)
assert box.npar == 0
system.add_particle(particle)
assert box.npar == 1
system.rm_particle(0)
assert box.npar == 0
system.add_particle(particle, box_id=0)
assert box.npar == 1
system.rm_particle(0, box_id=0)

# test the copy constructor
particle1 = mc.Particle("H", [0,0])
particle2 = mc.Particle(particle1)
box.add_particles([particle1, particle2])
assert box.npar == 2
assert system.ndim == 2
box.rm_particle(1)
box.rm_particle(0)
assert box.npar == 0

# add and remove single particle using box method
box.add_particle("H", [0,0,0])
assert box.npar == 1
assert system.ndim == 3
box.rm_particle(0)
assert box.npar == 0

# add and remove single particle using system method, default box
system.add_particle("H", [0,0,0])
assert box.npar == 1
system.rm_particle(0)
assert box.npar == 0

# add and remove single particle using system method, box_id used
system.add_particle("H", [0,0,0], box_id=0)
assert box.npar == 1
system.rm_particle(0, box_id=0)
assert box.npar == 0

# add list of particle objects
N = 10
particles = []
for i in range(N):
    particles.append(mc.Particle("H", [i, i, i]))
box.add_particles(particles)
assert box.npar == N
for i in range(N):
    box.rm_particle(0)
assert box.npar == 0
system.add_particles(particles)
assert box.npar == N
for i in range(N):
    system.rm_particle(0)
assert box.npar == 0
system.add_particles(particles, box_id=0)
assert box.npar == N
for i in range(N):
    system.rm_particle(0, box_id=0)
assert box.npar == 0

# add list of particle positions from fcc
n = 4
D = 3
N = (D+1) * n ** D
positions = mc.fcc(n, n * 1.7, D)
box.add_particles("H", positions)
assert box.npar == N
assert system.ndim == D
for i in range(N):
    box.rm_particle(0)
assert box.npar == 0
system.add_particles("H", positions)
assert box.npar == N
assert system.ndim == D
for i in range(N):
    system.rm_particle(0)
assert box.npar == 0
system.add_particles("H", positions, box_id=0)
assert box.npar == N
assert system.ndim == D
for i in range(N):
    system.rm_particle(0, box_id=0)
assert box.npar == 0

# read xyz-file
particles = mc.read_xyz("hydrogen.xyz")

# initialize from xyz-file
system.read_particles("hydrogen.xyz")

# test various dimensions
n = 4
N = 0
for D in [1, 2, 3]:
    N += (D+1) * n ** D
    positions = mc.fcc(n, n * 1.7, D)
    box.add_particles("H", positions)
    #assert box.npar == N
    assert system.ndim == D
