def sim_setup(bin):
    '''
    Conveniently setup simulation
    Earth radius
    returns:
     - sim (simulation pointer) 
     - ps (pointer to particle list)
     - rebx (pointer to reboundx link)
    '''
    import rebound
    import reboundx
    from src.utilities.constants import Rearth

    sim  = rebound.Simulation(bin)
    ps = sim.particles
    sim.integrator = "ias15" # choose ias15 integrator
    sim.testparticle_type = 1 # set masless particles to be non-interactive
    sim.dt = 1e4 # initial timestep, will be adjusted by integrator
    sim.N_active = sim.N # define active particles to be those in binary file

    ps[3].r = Rearth # set Earth Radius

    rebx = reboundx.Extras(sim)  # link simulation instance to reboundx
    rf = rebx.load_force("radiation_forces")  # load radiation forces package
    rebx.add_force(rf)  # include force in sim
    rf.params["c"] = 3.e8  # Speed of light in m/s
    ps["sun"].params["radiation_source"] = 1  # Assiging Sun as radiation source

    return sim, ps, rebx