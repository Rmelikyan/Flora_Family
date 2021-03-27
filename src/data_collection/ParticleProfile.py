from numpy.linalg import norm
class ParticleProfile:
    def __init__(self, p_hash, date, velocity, size, beta, is_collided, mean_anomaly):
        self.p_hash = p_hash
        self.date =  date
        self.vel = list(velocity)
        self.speed = round(norm(velocity), 2)
        self.size = size
        self.beta = beta
        self.is_collided = is_collided
        self.initial_mean_anomaly = mean_anomaly
    
    def collision(self, currTime):
        self.is_collided = True
        self.collision_date = currTime
