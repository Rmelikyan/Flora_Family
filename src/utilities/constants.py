import numpy as np

rad2Deg = 180/np.pi
au2m = 149597870700 # m/au
m2au = au2m**-1 # au/m
sec2year = 3.153600000e+7
secPerDay = 86400
G = 6.67408e-11

Msun = 1.3271244004193938e20/G
ss_masses = {
    'Sun':Msun,
    'Mercury':Msun/6023600,
    'Venus':Msun/408523.71, 
    "Earth":5.97237e24, 
    "Moon":7.346029776e22,
    "EMB": 5.97237e24 + 7.346029776e22,
    "Mars":Msun/3098708, 
    'Jupiter':Msun/1047.3486, 
    'Saturn':Msun/3497.898, 
    'Uranus':Msun/22902.98, 
    'Neptune':Msun/19412.24, 
    'Pluto':Msun/1.35e8,
    'Ceres':9.393e20, 
    'Pallas':2.05e20, 
    'Vesta':2.69e20
}

Rearth = 6.371e6  # m