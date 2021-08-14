from src.utilities.constants import rad2Deg
from src.utilities.misc import mod360

def invariable_inc(sim, p):
    from numpy import cross, dot, arccos
    from numpy.linalg import norm
    sim_ang_mom = sim.calculate_angular_momentum()
    p_ang_mom = cross(p.xyz, p.vxyz)

    inc = arccos(dot(sim_ang_mom, p_ang_mom)/(norm(sim_ang_mom)*norm(p_ang_mom)))
    return inc*rad2Deg

def equinoctal_2_eccliptic(a, g=None, f=None, k=None, h=None, L=None):
    import numpy as np
    if len(a) == 6:
        L = a[5]
        h = a[4]
        k = a[3]
        f = a[2]
        g = a[1]
        a = a[0]
    e = np.sqrt(f**2 + g**2)
    i = mod360(np.arctan2(2*np.sqrt(h**2 + k**2), 1 - h**2 - k**2) *rad2Deg)
    om = mod360(np.arctan2(g*h - f*k, f*h + g*k) * rad2Deg)
    OM = mod360(np.arctan2(k,h) * rad2Deg)
    true_anomaly = mod360((L - (OM + om)))
    return [a, e, i, om, OM, true_anomaly]

