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
    'Neptune':Msun/19412.24
}

small_body_masses = {
    'Pluto':Msun/1.35e8,
    'Ceres':9.38e20, 
    'Pallas':2.11e20, 
    'Vesta':2.59e20,
    'A847 NA': 1.34e19, # Hebe
    'A847 PA': 1.79e19, # Iris
    'A848 HA': 1.16e19, # Metis
    'Parthenope': 5.86e18,
    'Victoria': 3.56e18,
    'Melpomene': 3.00e18,
    'Fortuna': 8.35e18,
    'Massalia': 5.67e18,
    'Lutetia': 2.06e18,
    'Euterpe': 4.35e18,
    'Harmonia': 1.99e18,
    'Isis': 1.86e18,
    'Nemausa': 1.99e18,
    '105': 2.89e18, # Artemis name gets confused with spacecraft
    'Athamantis': 1.89e18,
    'Zelinda': 1.35e18 
}

Rearth = 6.371e6  # m