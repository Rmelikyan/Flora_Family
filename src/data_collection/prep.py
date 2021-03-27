from src.utilities.misc import zero_to_2pi
from src.utilities.constants import rad2Deg, m2au
import numpy as np
from scipy.stats import  circstd

def sigma_m(orbits):
    if len(orbits) < 1:
        return 0
    big_ms = [zero_to_2pi(o.M) for o in orbits]
    std_ms = circstd(big_ms)*rad2Deg
    return std_ms

def dsh_of_stream(ps):
    if len(ps) < 2:
        return 0
    parent_body_repr = ps[-1]
    a1 = parent_body_repr.a * m2au
    e1 = parent_body_repr.e
    inc1 = parent_body_repr.inc
    pomega1 = zero_to_2pi(parent_body_repr.pomega)
    Omega1 = zero_to_2pi(parent_body_repr.Omega)

    dshs = []
    for p in ps[:-1]:
        a2 = p.a * m2au
        e2 = p.e
        inc2 = p.inc
        pomega2 = zero_to_2pi(p.pomega)
        Omega2 = zero_to_2pi(p.Omega)

        dsh = D_sh(a1,a2,e1,e2,inc1,inc2,pomega1,pomega2,Omega1,Omega2)
        dshs.append(dsh)
    return [np.mean(dshs), np.std(dshs)]


def D(a1,a2,e1,e2,inc1,inc2):
    ''' Calculate D parameter'''
    q1 = a1*(1-e1) #Periapsis
    q2 = a2*(1-e2)
    
    D2 = (e2 - e1)**2 + (q2 - q1)**2 + (2*np.sin((inc2-inc1)/2))**2 #D-Criteron
    return np.sqrt(D2)

def D_sh(a1,a2,e1,e2,inc1,inc2,pomega1,pomega2,Omega1,Omega2):
    ''' Calculate Dsh parameter'''
    Dc = D(a1,a2,e1,e2,inc1,inc2) #D-Criteron
    D2 = Dc**2
    
    D_sh2 = D2 + np.sin(inc1)*np.sin(inc2)*(2*np.sin((Omega2-Omega1)/2))**2 + ((e1+e2)/2 * 2*np.sin((pomega2-pomega1)/2))**2 #D_sh-Criteron
    return np.sqrt(D_sh2)

def helioStateVector(p, sim):
    ''' Returns helocentric state vector in au and au/year'''
    sim.move_to_hel()
    return p.xyz + p.vxyz

def stream(ps):
    if len(ps) < 1:
        return None
    elements = [[], [], [], [], []]
    for p in ps:
        elements[0].append(p.a)
        elements[1].append(p.e)
        elements[2].append(p.inc)
        elements[3].append(zero_to_2pi(p.omega))
        elements[4].append(zero_to_2pi(p.Omega))

    return [[np.mean(element), np.std(element)] for element in elements]


def getSolarSystemDict(solar_system, sim):
    names = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto", "Ceres", "Pallas", "Vesta"]
    ss_dict = {}
    for i, planet in enumerate(solar_system):
        name = names[i]
        ss_dict[name] = helioStateVector(planet, sim)
    return ss_dict