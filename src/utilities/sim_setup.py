from numpy import log
from src.utilities.constants import au2m, rad2Deg

def sim_setup(bin, log_msgs = None):
    '''
    Conveniently setup simulation
    returns:
     - sim (simulation pointer)
    '''
    import rebound
    if log_msgs is None:
        log_msgs = []
    sim  = rebound.Simulation(bin)
    ps = sim.particles
    sim.integrator = "mercurius" # choose mercurius integrator
    log_msgs.append("sim.integrator = 'mercurius'\n")
    sim.testparticle_type = 1 # set masless particles to be non-interactive
    log_msgs.append("sim.testparticle_type = 1\n")
    sim.dt = ps[1].P* 0.05 # initial timestep, will be adjusted by integrator
    log_msgs.append("sim.dt = {} (s)\n".format(sim.dt))
    sim.N_active = sim.N # define active particles to be those in binary file

    a = 329320844896.02423 # m 
    e = 0.15585014794810317
    i = 5.889091687694949/rad2Deg
    Om = 110.87633970146754/rad2Deg
    om = 285.50181319194064/rad2Deg
    M  = 15.638508490829182/rad2Deg
    sim.add(a = a, e = e, inc = i, Omega = Om, omega = om, M = M, hash = 'Flora')
    log_msgs.append("Flora added at pos (m):" +
                    "{}, vel (m/s): {}\n".format(ps[-1].xyz, ps[-1].vxyz))
    return sim