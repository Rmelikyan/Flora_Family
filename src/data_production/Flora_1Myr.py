import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from datetime import datetime
import rebound
import numpy
from src.utilities.sim_setup import sim_setup
from src.utilities.constants import sec2year

start_time = datetime.now()
log_msgs = ["Running {} at {}\n".format(__name__, start_time)]

sim_path = 'data/sim_inits/2020_12_17_sim.bin'
log_msgs.append("Sim initialized by sim_setup.py using init_file: {}\n".format(sim_path))
sim = sim_setup(sim_path, log_msgs)

file_name = 'Flora_1Myr_mercurius.bin'
archive_path = 'data/sim_archives'
save_name = archive_path+'/'+file_name
tmax = 1e6*sec2year
log_msgs.append('tmin = {} --> tmax = {} yrs\n'.format(sim.t, tmax/sec2year))

archive_interval = 50*sec2year

sim.automateSimulationArchive(save_name, archive_interval, deletefile=True)
log_msgs.append('Archival Interval: {}\n'.format(archive_interval))
log_msgs.append('Archived File Path: {}\n'.format(save_name))

print('Starting!')
sim.integrate(tmax)

end_time = datetime.now()
log_msgs.append('End Time: {}\n'.format(end_time))

print(end_time - start_time)
log_msgs.append('Total Run Time: {}\n'.format(end_time - start_time))

with open(save_name[:-3] + 'log', 'w') as log:
    for msg in log_msgs:
        log.write(msg)