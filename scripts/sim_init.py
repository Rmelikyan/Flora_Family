import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import rebound
from datetime import datetime
from src.utilities.constants import ss_masses, small_body_masses

def sim_init(date):
    log_msgs = ["Initializing Sim at date: {}\n".format(date)]
    formatted_date = str(date)[:10].replace('-', '_')

    sim = rebound.Simulation()
    sim.units = ('s', 'm', 'Kg')
    log_msgs.append("Sim units used: (s, m, Kg)\n")

    perturbers = ["Sun", "Mercury", "Venus", "Earth", "Mars",
                  "Jupiter", "Saturn", "Uranus", "Neptune"]
    small_body_names = list(small_body_masses.keys())

    perturbers += small_body_names
    log_msgs.append("Perturbers included: {}\n".format(str(perturbers)))


    for name in perturbers:
        if name == "Earth" and "Moon" not in perturbers:
            mass_key = "EMB"
        else:
            mass_key = name

        sim.add(name, date = date)
        sim.particles[-1].hash = name
        if mass_key in ss_masses.keys():
            sim.particles[-1].m = ss_masses[mass_key]
        elif mass_key in small_body_masses.keys():
            sim.particles[-1].m = small_body_masses[mass_key]
    log_msgs.append("Masses as defined in constants.py used\n")
    log_msgs.append("NOTE: If Earth included but not Moon then EMB mass used\n")

    sim.move_to_hel()
    log_msgs.append("Bodies have been shifted to heliocentric positions\n")


    save_name='data/sim_inits/{}_sim_SMALLBODIES.bin'.format(formatted_date)
    sim.save(save_name)

    log_msgs.append("sim.bin file located at: {}\n".format(save_name))
    log_msgs.append("Completed: {}\n".format(datetime.now()))

    with open(save_name[:-3] + 'log', 'w') as log:
        for msg in log_msgs:
            log.write(msg)

    return

if __name__ == "__main__":
    date = datetime(2020, 12, 17, 0, 0, 0)
    sim_init(date)