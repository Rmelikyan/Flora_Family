def chunk (f_path, production_years, test_name, ppw = 1, notebook = False):
    '''
    Runs simulation and saves data to f_path
    production_years defines range of years for particle production
    '''
    import json
    from datetime import datetime, timedelta

    import numpy as np
    import rebound
    import reboundx
    from src.data_collection.my_converter import my_converter
    from src.data_collection.prep import (dsh_of_stream, helioStateVector,
                                          sigma_m, stream)
    from src.utilities.constants import Rearth, m2au
    from src.utilities.ephem_reader import read_jpl_ephemeris
    from src.utilities.misc import T_print, normDif
    from src.utilities.particle_production import fullEjection
    from src.utilities.sim_setup import sim_setup
    begin = datetime.now()

    def myResolution(a, b):
        ps = sim.particles
        swap = 0 
        i = b.p1
        j = b.p2
        if i > j:
            swap = 1
            i = b.p2
            j = b.p1

        collTime = initTime + timedelta(0, sim.t)
        dist = normDif(ps[i].xyz, ps[j].xyz)

        if ps[i].m + ps[j].m == 0:
            return 0
        elif dist > Rearth:
            if ps[j].hash.value in data[collTime.year]["Particles"].keys():
                prev_dist = data[collTime.year]["Particles"][ps[j].hash.value]["Dist"]
                if dist < prev_dist:
                    data[collTime.year]["Particles"][ps[j].hash.value] = {"State" : helioStateVector(ps[j], sim), "Time" : collTime, "Dist" : dist}
            else:
                data[collTime.year]["Particles"][ps[j].hash.value] = {"State" : helioStateVector(ps[j], sim), "Time" : collTime, "Dist" : dist}
            return 0
        
            
        ps[i].x = (ps[i].x * ps[i].m + ps[j].x * ps[j].m)/(ps[i].m + ps[j].m)
        ps[i].y = (ps[i].y * ps[i].m + ps[j].y * ps[j].m)/(ps[i].m + ps[j].m)
        ps[i].z = (ps[i].z * ps[i].m + ps[j].z * ps[j].m)/(ps[i].m + ps[j].m)
        ps[i].vx = (ps[i].vx * ps[i].m + ps[j].vx * ps[j].m)/(ps[i].m + ps[j].m)
        ps[i].vy = (ps[i].vy * ps[i].m + ps[j].vy * ps[j].m)/(ps[i].m + ps[j].m)
        ps[i].vz = (ps[i].vz * ps[i].m + ps[j].vz * ps[j].m)/(ps[i].m + ps[j].m)
        ps[i].m = ps[i].m + ps[j].m
        ps[i].r  = pow(pow(ps[i].r,3.)+pow(ps[j].r,3.),1./3.)
        
        particle_library[ps[j].hash.value].collision(collTime)

        return 1 if swap else 2
    if notebook:
        sim, ps, rebx = sim_setup('../data/sim_inits/1900_01_01_sim.bin')
    else:
        sim, ps, rebx = sim_setup('data/sim_inits/1900_01_01_sim.bin')
    sim.N_active = sim.N - 1
    ps[3].r *= 500

    sim.collision = "direct"
    sim.collision_resolve = myResolution

    startYear = 1900
    endYear = 2029

    if production_years[0] > endYear:
        return 0

    if notebook:
        Time, X, Y, Z, VX, VY, VZ = read_jpl_ephemeris('../data/ephemerides/1900_01_01Apophis.txt')
    else:
        Time, X, Y, Z, VX, VY, VZ = read_jpl_ephemeris('data/ephemerides/1900_01_01Apophis.txt')
    

    # ------------------ Path Support ------------------
    f_path.mkdir(parents=True, exist_ok=True)
    # -----------------------------------------------------------
    
    initTime = Time[0]
    startTime = Time[0]
    for t in Time:
        if t.year == startYear:
            startTime = t
            break
    endTime = datetime(endYear, 4, 30, 0, 0, 0)

    GT_particledays_1ppw = 159397882
    particle_days = 0
    
    data = {}
    particle_library = {}
    collision_messages = []
    data_count = 0
    prev_particle_num = sim.N
    num_collisions = 0
    num_ps = 0

    data[startYear] = {"Particles": {}}
    data[startYear]["Dsh"] = dsh_of_stream(ps[13:])
    data[startYear]["Stream"] = stream(ps[13:])
    data[startYear]["Sigma_M"] = sigma_m(sim.calculate_orbits(primary=ps[0])[13:])

    currTime = startTime

    while currTime <= endTime:

        # ----------------------- Integration -----------------------
        seconds = (currTime-initTime).total_seconds()
        sim.integrate(seconds, exact_finish_time=1)
        sim.move_to_hel()
        # -----------------------------------------------------------

        # ------------------ Data Collection logic ------------------
        if currTime.month == 1 and currTime.day == 1:
            key = currTime.year
            data[key] = {"Particles": {}}
            data[key]["Dsh"] = dsh_of_stream(ps[13:])
            data[key]["Stream"] = stream(ps[13:])
            data[key]["Sigma_M"] = sigma_m(sim.calculate_orbits(primary=ps[0])[13:])
        if currTime.month == 4  and currTime.day == 30:
            data_count = sum([len(data[key]["Particles"]) for key in data.keys()])
        # -----------------------------------------------------------

        # -------------------- Collision Capture --------------------
        if (prev_particle_num - sim.N) > 0:
            curr_collisions = (prev_particle_num - sim.N)
            
            coll_msg = "\n{} Collisions on {}".format(curr_collisions, currTime)
            collision_messages.append(coll_msg)

            num_collisions += curr_collisions

        prev_particle_num = sim.N 
        # -----------------------------------------------------------

        # -------------------- Particle Creation --------------------
        if currTime.year in production_years and currTime in Time:
            i = Time.index(currTime)
            state = [X[i], Y[i], Z[i], VX[i], VY[i], VZ[i]]
            p_profiles = fullEjection(ppw, state, sim, rebx, currTime)
            for profile in p_profiles:
                particle_library[profile.p_hash.value] = profile
                num_ps += 1
        # -----------------------------------------------------------

        if notebook:
            print(currTime, '\t', num_ps, '\t', data_count, '\t', num_collisions, '\t', data[currTime.year]["Sigma_M"], end='\r')

        particle_days += (num_ps) 
        currTime += timedelta(1) # Next day

    # ----------------------- Timing Analysis -----------------------
    dt = datetime.now() - begin # total integration time
    p_days_per_sec = particle_days/dt.total_seconds() # integration velocity
    expected_gt_comp_time = timedelta(0, GT_particledays_1ppw * ppw/p_days_per_sec)
    # ---------------------------------------------------------------
    
    # ------------------------ File Logging -------------------------
    with open(f_path / '{}.txt'.format(test_name), "w") as txt:

        T_print("Total Time Range: {} - {}".format(startTime, endTime), txt)

        T_print("\nProduction Range: {} - {}".format(production_years[0], production_years[-1]), txt)

        T_print("\nParticle sizes range from .1-2 (cm) diameter according to LogNormal Distribution", txt)

        T_print("\nParticle ejection velocities range from .2-1 (m/s) and are distributed uniformly", txt)

        message = "\nEnd Date: {}\tTotal Particles: {}\t Close Encounters {}\t Internal Collisions {}".format(currTime, num_ps, data_count, num_collisions)
        T_print(message, txt)
        
        T_print("\nTotal integration time (H:M:S): {}".format(dt), txt)

        T_print('\nParticle Days Simulated: {}'.format(particle_days), txt)

        T_print("ParticleDays per Second: {}".format(p_days_per_sec), txt)

        T_print("Expected total computational time for 1 core: {}".format(expected_gt_comp_time), txt)

        T_print("\n----Collision Messages----", txt)
        for msg in collision_messages:
            T_print(msg, txt)

    jsonOut = f_path / "{}.json".format(test_name)
    with open(jsonOut, "w") as fp:
        json.dump(data, fp, indent=2, default=my_converter)  # Dump dict to jsonOut file

    jsonOut = f_path / "{}_particles.json".format(test_name)
    with open(jsonOut, "w") as fp:
        json.dump(particle_library, fp, indent=2, default=my_converter)  # Dump dict to jsonOut file
    # ---------------------------------------------------------------

    if notebook:
        return [particle_days, dt, num_ps, num_collisions, sim]
    else:
        return [particle_days, dt, num_ps, num_collisions]
