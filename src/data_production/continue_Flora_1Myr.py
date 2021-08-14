import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from datetime import date, datetime
import rebound
import numpy
from src.utilities.sim_setup import sim_setup
from src.utilities.constants import sec2year

start_time = datetime.now()

sim_init_path = 'data/sim_inits/2020_12_17_sim.bin'

file_name = 'Flora_1Myr.bin'
archive_path = 'data/sim_archives'

sim = rebound.Simulation(archive_path+'/'+file_name)
sim.integrator = 'mercurius'
sim.dt = sim.particles[1].P* 0.05

archive_interval = 50*sec2year
tmax = 1e6*sec2year
sim.automateSimulationArchive(archive_path+'/'+file_name, archive_interval, deletefile=False)
print('Starting!')
while sim.t < tmax:
    t_start = datetime.now()
    sim.integrate(sim.t+ 10000*sec2year)
    t_end = datetime.now()
    print("Sim at t = {} years, last thousand years took: {} (s)".format(sim.t/sec2year, t_end-t_start), end = '\r')

end_time = datetime.now()

print(end_time - start_time)