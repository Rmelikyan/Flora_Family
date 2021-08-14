import rebound
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from src.utilities.particle_production import EjectionVector


sim = rebound.Simulation('/Users/bethclark/Projects/Flora_Family/data/sim_inits/Flora_Init_1.bin')
sim.remove(hash='Mercury') # this guy slows things down
sim.move_to_hel()
ps = sim.particles
eject_vectors = [EjectionVector(100) for i in range(20)] # initial V vectors for 20 ps
x, y, z = ps["Flora"].xyz
vx, vy, vz = ps["Flora"].vxyz
sim.remove(hash="Flora") # this guy is unecessary
for V in eject_vectors:
    sim.add(x=x, y=y, z=z,vx=vx+V[0], vy=vy+V[1], vz=vz+V[2], r = 5000)

sim2 = rebound.Simulation()
sim2.units = ("m", "kg", 's')
sim2.integrator = "WHFast"
sim2.dt = 1e6
for p in ps:
    sim2.add(p)
#     x, y, z = p.xyz
#     vx, vy, vz = p.vxyz
#     m = p.m
#     r = p.r
#     sim2.add(x=x, y=y, z=z, vx=x, vy=y, vz=z, m=m, r=r)

del sim
del sim2
print('done')