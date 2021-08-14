/**
 * Radiation forces
 *
 * This example provides an implementation of the 
 * Poynting-Robertson effect. The code is using the IAS15 integrator
 * which is ideally suited for this velocity dependent force.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void yarko_da(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);

double sec2year = 31536000.0;
double rad2Deg = 57.29577951308232;
double years = 1e4;
double dadts[20] = {2e-5, -2e-5, 2e-5, -2e-5, 2e-5, -2e-5, 2e-5, -2e-5, 2e-5, -2e-5, 2e-5, -2e-5, 2e-5, -2e-5, 2e-5, -2e-5, 2e-5, -2e-5, 2e-5, -2e-5};
double au2m = 149597870700;

int main(int argc, char* argv[]){
    if(argc != 2){
        printf("Must include 'path_to_sim' as command line inputs!\nProgram Terminated\n");
        exit(0);
    }

    double tmax = (years + 1) * sec2year;

    char* path_to_sim = argv[1];

    struct reb_simulationarchive* sa = reb_open_simulationarchive(path_to_sim);
    if (sa==NULL){
        printf("Can not open sim file.\nProgram Temrinating\n");
        exit(0);
    }
    // Get a simulation from the file (if possible, otherwise NULL is returned)
    struct reb_simulation* r = reb_create_simulation_from_simulationarchive(sa,0);
    // Whenever you've opened a SimulationArchive and don't need it anymore, close it.
    reb_close_simulationarchive(sa);
    // Check if we were successful

    remove(path_to_sim); // remove init file. To Be Replaced with archive

    reb_simulationarchive_automate_interval(r,path_to_sim,500*sec2year);
    r->force_is_velocity_dependent  = 1;
    r->additional_forces            = yarko_da;    // setup callback function for velocity dependent forces
    r->heartbeat                    = heartbeat;
    r->N_active                     = 8;
    reb_integrate(r, tmax);

    // reb_free_pointers(r);
    // reb_free_simulation(r);
    exit(0);
}

void yarko_da(struct reb_simulation* r){
    struct reb_particle* particles = r->particles;
    const int N = r->N;
    const struct reb_particle star = particles[0];                // cache
    int j = 0;
#pragma omp parallel for
    for (int i=0;i<N;i++){
        const struct reb_particle p = particles[i];             // cache
        if (p.m!=0.) continue;                         // only dust particles feel radiation forces
        const double prx  = p.x-star.x;
        const double pry  = p.y-star.y;
        const double prz  = p.z-star.z;
        const double pr   = sqrt(prx*prx + pry*pry + prz*prz);         // distance relative to star
        
        const double prvx = p.vx-star.vx;
        const double prvy = p.vy-star.vy;
        const double prvz = p.vz-star.vz;
        const double v2   = prvx*prvx + prvy*prvy + prvz*prvz;

        const double gm         = r->G*star.m;
        const double energy     = 0.5*v2 - gm/pr;
        const double a          = -0.5*gm/energy;
        const double auPmyr2mPs = au2m/sec2year/1e6;
        const double dadt       = dadts[j]*auPmyr2mPs;
        const double k          = 0.5*dadt*gm/(a*a);

        particles[i].ax += k * prvx/v2;
        particles[i].ay += k * prvy/v2;
        particles[i].az += k * prvz/v2;
        j+=1;
    }
}

void heartbeat(struct reb_simulation* r){
    if(reb_output_check(r, 1e3*sec2year)){                        // print some information to screen
        reb_output_timing(r, years * sec2year);
    }
}
