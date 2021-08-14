from ctypes import ArgumentError
import subprocess
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

def run_10km_yarko(run_num):
    import rebound
    from pathlib import Path
    from src.utilities.particle_production import EjectionVector

    # -------------- Set Up Sim Init File --------------
    sim = rebound.Simulation('/Users/bethclark/Projects/Flora_Family/data/sim_inits/Flora_Init_1.bin')
    sim.remove(hash='Mercury') # this guy slows things down
    sim.move_to_hel()
    ps = sim.particles
    eject_vectors = [EjectionVector(100) for i in range(20)] # initial V vectors for 20 ps
    x, y, z = ps["Flora"].xyz
    vx, vy, vz = ps["Flora"].vxyz
    sim.remove(hash="Flora") # this guy is unecessary
    for V in eject_vectors:
        sim.add(x=x, y=y, z=z,vx=vx+V[0], vy=vy+V[1], vz=vz+V[2], r = 5000)

    sim2 = rebound.Simulation()
    sim2.units = ("m", "kg", 's')
    sim2.integrator = "WHFast"
    sim2.dt = 1e6
    for p in ps:
        sim2.add(p)

    init_path = Path('/Users/bethclark/Projects/Flora_Family/data/sim_archives')
    test_name = 'Flora_10km_run_' + run_num

    init_file_path = (init_path / test_name).__str__()

    sim2.save(init_file_path)
    
    del sim
    del sim2

    program_path = './src_C/10km_yarko/rebound'

    subprocess.run([program_path, init_file_path])
    print('\nfinished C... now returning')
    return

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise ArgumentError

    run_num = sys.argv[1]
    run_10km_yarko(run_num)