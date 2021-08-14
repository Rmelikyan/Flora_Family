import numpy as np
def pm(a, b):
    return (a-b, a+b)

def zero_to_2pi(angle):
    return angle if angle >= 0 else angle + 2*np.pi

def T_print(message, file):
    """ prints and writes message to file"""
    print(message)
    file.write(str(message) + '\n')

def normDif(a1, a2):
    norm = np.linalg.norm
    if len(a1) is not len(a2):
        raise Exception("Bad Arrays")

    return(norm([a1[i]-a2[i] for i in range(len(a1))]))

def mod360(ang):
    """ Transforms degree angle of arbitrary range to 0 - 360 """
    while ang < 0:
        ang+=360
    return ang%360
        