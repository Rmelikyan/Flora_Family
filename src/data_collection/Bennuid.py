import numpy as np
from numpy.linalg import norm

class Bennuid:
    ''' A convience class for Bennuid Data analysis
    
        A mild copy of the rebound particle instance without Sim dependance
    '''
    au2m = 149597870700 # m/au
    sec2year = 3.153600000e+7
    Msun = 1.988475415966536e+30
    Gsi = 6.67408e-11
    
    def __init__(self, data, units, refFrame, date):
        self.x  = data[0]
        self.y  = data[1]
        self.z  = data[2]
        self.vx = data[3]
        self.vy = data[4]
        self.vz = data[5]
        self.a = None
        self.e = None
        self.inc = None
        self.omega = None
        self.Omega = None
        self.units = units
        self.refFrame = refFrame
        self.date = date
    
    def convertHelioCentric(self, earth):
        if (self.units != earth.units):
            raise Exception("Particle and Earth not of same Units")
        if (self.date != earth.date):
            raise Exception("Particle and Earth not of same Date")
        if (self.refFrame == 'Helio'):
            raise Exception("Particle is already in Heliocentric Coordinates")
        if (earth.refFrame != 'Helio'):
            raise Exception("Earths must be in Helio")
        self.x  += earth.x
        self.y  += earth.y
        self.z  += earth.z
        self.vx += earth.vx
        self.vy += earth.vy
        self.vz += earth.vz
        self.refFrame = "Helio"
        
    def convertGeoCentric(self, earth):
        if (self.units != earth.units):
            raise Exception("Particle and Earth not of same Units")
        if (self.date != earth.date):
            raise Exception("Particle and Earth not of same Date")
        if (self.refFrame == 'Geo'):
            raise Exception("Particle is already in Geocentric Coordinates")
        if (earth.refFrame != 'Helio'):
            raise Exception("Earths must be in Helio")
        self.x  -= earth.x
        self.y  -= earth.y
        self.z  -= earth.z
        self.vx -= earth.vx
        self.vy -= earth.vy
        self.vz -= earth.vz
        self.refFrame = "Geo"
        
    def calcElements(self, M = Msun, G = Gsi):
        if (G == self.Gsi and self.units != 'SI'):
            raise Exception("Unit Error")
        mu = M * G
        hVector = np.cross(self.xyz(), self.vxyz())
        self.h = hVector
        eVector = np.cross(self.vxyz(), hVector/mu) - self.xyz()/norm(self.xyz())
        nVector = np.cross([0,0,1], hVector)
        self.inc = np.arccos(hVector[2]/norm(hVector))
        self.e = norm(eVector)
        if (nVector[1] >= 0):
            self.Omega = np.arccos(nVector[0]/norm(nVector))
        else:
            self.Omega = 2*np.pi - np.arccos(nVector[0]/norm(nVector))
        if (eVector[2] >= 0):
            self.omega = np.arccos(np.dot(nVector,eVector)
                                    /(norm(nVector)*norm(eVector)))
        else:
            self.omega = 2*np.pi - np.arccos(np.dot(nVector,eVector)
                                                /(norm(nVector)*norm(eVector)))
        if (np.dot(self.xyz(),self.vxyz()) >= 0):
            self.f = np.arccos(np.dot(eVector,self.xyz())/(self.e*norm(self.xyz())))
        else:
            self.f = 2*np.pi - np.arccos(np.dot(eVector,self.xyz())/(self.e*norm(self.xyz())))
        
        self.a = (2/norm(self.xyz()) - norm(self.vxyz())**2/mu)**-1
        if self.a > 0:
            self.P = 2*np.pi*np.sqrt(self.a**3/mu)
        else:
            self.P = -1
        self.n = 2*np.pi/self.P

        
    def xyz(self):
        return [self.x, self.y, self.z]
    
    def vxyz(self):
        return [self.vx, self.vy, self.vz]

    def convertSI(self):
        self.x *= Bennuid.au2m
        self.y *= Bennuid.au2m
        self.z *= Bennuid.au2m
        self.vx *= Bennuid.au2m/Bennuid.sec2year
        self.vy *= Bennuid.au2m/Bennuid.sec2year
        self.vz *= Bennuid.au2m/Bennuid.sec2year
        self.units = 'SI'

    def convertAU(self):
        self.x /= Bennuid.au2m
        self.y /= Bennuid.au2m
        self.z /= Bennuid.au2m
        self.vx /= Bennuid.au2m/Bennuid.sec2year
        self.vy /= Bennuid.au2m/Bennuid.sec2year
        self.vz /= Bennuid.au2m/Bennuid.sec2year
        self.units = 'AU'

    def convertG1(self):
        self.x /= Bennuid.au2m
        self.y /= Bennuid.au2m
        self.z /= Bennuid.au2m
        self.vx /= Bennuid.au2m/Bennuid.sec2year*2*np.pi
        self.vy /= Bennuid.au2m/Bennuid.sec2year*2*np.pi
        self.vz /= Bennuid.au2m/Bennuid.sec2year*2*np.pi
        self.units = 'G1'

    def BplaneAnalysis(self, Mearth, Rearth):
        if (self.refFrame == 'Helio'):
            raise Exception("Particle must be in Geocentric Frame")
        if self.units != 'SI':
            self.convertSI()
        if self.a is None:
            self.calcElements(M=Mearth)
        mu = self.Gsi * Mearth
        V_inf = np.sqrt(-mu/self.a)
        B = np.linalg.norm(self.h)/V_inf
        CaptureDistance = Rearth*np.sqrt(1+2*mu/(Rearth*V_inf**2))
        return B, CaptureDistance, bool(B < CaptureDistance)

    def Bvector(self, Mearth):
        if (self.refFrame == 'Helio'):
            raise Exception("Particle must be in Geocentric Frame")
        if self.units != 'SI':
            self.convertSI()
        if self.a is None:
            self.calcElements(M=Mearth)
        mu = self.Gsi * Mearth
        V_inf = np.sqrt(-mu/self.a)
        e_vec = np.cross(self.vxyz(),self.h)/mu - self.xyz()/norm(self.xyz())
        P_hat = e_vec/self.e
        Q_hat = np.cross(self.h,P_hat)/norm(self.h)
        S_hat = 1/self.e*P_hat + np.sqrt(self.e**2-1)/self.e*Q_hat
        B_vec = np.cross(S_hat, self.h)/V_inf
        return B_vec, S_hat

    def OpikAnalysis(self, Vpl, Mearth):
        ''' Calculate the Moid xi and timing zeta
        
            Vpl should be the heliocentric velocity vector of Earth
        '''
        if (self.refFrame == 'Helio'):
            raise Exception("Particle must be in Geocentric Frame")
        if self.a is None:
            self.calcElements(M=Mearth)
        if self.units != 'SI':
            self.convertSI()
        B_vec, S_hat = self.Bvector(Mearth)
        XI = np.cross(Vpl, S_hat)
        XI = XI/norm(XI)
        ZETA = np.cross(-S_hat, XI)
        xi = np.dot(B_vec, XI)
        zeta = np.dot(B_vec, ZETA)
        if (abs(xi) > norm(B_vec)):
            print("Why is moid greater than impact estimate?...")
        return xi, zeta
        
