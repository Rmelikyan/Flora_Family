import random
import numpy as np
from scipy.stats import lognorm
import rebound
import datetime
from src.data_collection.ParticleProfile import ParticleProfile

def randVelocity():
    ''' get random speed between the range of likely ejection speeds (.2-1 m/s)'''
    minVelocity = .2 #Escape velocity defines lower limit: .2 m/s
    maxVelocity = 1 #Fastest observed particles: 1 m/s
    V = random.randrange(int(1e4*minVelocity), int(1e4*maxVelocity))/1e4 # Randomly chosen velocity
    return V

def EjectionVector(velocity):
    ''' Returns random vectorized components for a given ejection speed'''
    vector = np.random.rand(3)-.5 # random cartesion vector of positive and negative numbers
    vector = vector/np.linalg.norm(vector)
    v_vector = vector*velocity
    return list(v_vector)

def randRadius():
    '''Using Bennu particle production distribution to generate particle sizes'''
    lognorm_dist = lognorm(0.5288479765345986, -0.7695466875271932, 3.351233773990743)
    hist = np.histogram(lognorm_dist.rvs(size = 1), bins=np.linspace(0,10, 26), range=None, normed=None, weights=None, density=None)
    try:
        size_not_scaled = round(hist[1][hist[0].nonzero()[0]][0], 3)
    except:
        size_not_scaled = 2.6
    size = size_not_scaled * 1/4 * 1/100 + .001 # scaled
    return round(size/2, 5)

def calcBeta(grain_radius, sim, rebx):
    ''' Calculates and applies beta parameter to last created partilce 
    
        grain_radius (float): typically provided by randRadius'''
    density = 2000. #kg/m^3 = 2g/cc !!! Assumption
    Q_pr = 1. #Reflection effeciency
    luminosity = 3.85e26 #Watts
    Msun = sim.particles[0].m
    beta = rebx.rad_calc_beta(sim.G, 3.e8, Msun, luminosity, grain_radius, density, Q_pr) #Rebx function
    sim.particles[-1].params["beta"] = beta #Assign beta to most recent particle in simulation
    return beta

def fullEjection(numEjected, state, sim, rebx, curr_date, speed = None, radius = None):
    ''' Particle Ejection function which applies SRPR AND Ejection Velocity'''

    p_profiles = []
    rad2Deg = 180/np.pi

    for _ in range(numEjected):
        if (speed == None):
            v_vector = EjectionVector(randVelocity())

        else:
            v_vector = EjectionVector(speed)

        p_hash = rebound.hash(datetime.datetime.now().__str__() + str(v_vector)) # assign random hash

        if (radius == None):
            radius = randRadius()

        x, y, z, vx, vy, vz = state # unpack particle state vector

        sim.add(x=x, y=y, z=z,vx=vx+v_vector[0], vy=vy+v_vector[1], vz=vz+v_vector[2], r = radius, hash = p_hash) #adds ejection vector
        
        beta = calcBeta(radius, sim, rebx)

        m = sim.particles[-1].M # record initial mean anomaly
        m = m * rad2Deg if m > 0 else (m + 2*np.pi) * rad2Deg # convert to 0-2pi

        p_profiles.append(ParticleProfile(p_hash, curr_date, v_vector, radius, beta, False, m))
        
    return p_profiles