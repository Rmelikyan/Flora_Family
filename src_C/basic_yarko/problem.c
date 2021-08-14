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

double betaparticles = 3.321012071791588e-09;     // beta parameter, defined as the ratio of radiation pressure over gravity
double sec2year = 31536000.0;
double rad2Deg = 57.29577951308232;
double years = 1e6;
double dadts[5] = {1, .5, .1, 0.05, 0.01};
double au2m = 149597870700;

int main(int argc, char* argv[]){
    double tmax = years * sec2year;
    char* init_file = "/Users/bethclark/Projects/Flora_Family/data/sim_inits/2020_12_17_sim.bin";
    char* filename = "Yarko_1Myr_C_WHF.bin";

    struct reb_simulationarchive* sa = reb_open_simulationarchive(NULL);
    if (sa==NULL){
        printf("Can not open file.\n");
    }
    // Get a simulation from the file (if possible, otherwise NULL is returned)
    struct reb_simulation* r = reb_create_simulation_from_simulationarchive(sa,-1);
    // Whenever you've opened a SimulationArchive and don't need it anymore, close it.
    reb_close_simulationarchive(sa);
    // Check if we were successful
    if (r==NULL){
        printf("No simulation archive found. Creating new simulation.\n");
        r = reb_create_simulation();
        // setup constants
        r->G                            = 6.67408e-11;
        r->integrator                   = REB_INTEGRATOR_WHFAST;
        r->dt                           = 1e6;            // initial timestep
        // r->ri_ias15.epsilon             = 1e-4;            // accuracy parameter
        r->N_active                     = 1;             // the star is the only massive particle
        r->force_is_velocity_dependent  = 1;
        r->additional_forces            = yarko_da;    // setup callback function for velocity dependent forces
        r->heartbeat                    = heartbeat;
        
        // star is at rest at origin
        struct reb_particle star = {0};
        star.m  = 1.3271244004193938e20/r->G;
        reb_add(r, star);

        // dust particles are initially on a circular orbit
        
        for (int i=0; i<5; i++){
            double offset = 1e3*rand()/RAND_MAX;
            struct reb_particle primary = r->particles[0];
            double m = 0.;
            double a = 329320844896.02423 + offset; // m 
            double e = 0.15585014794810317;
            double inc = 5.889091687694949/rad2Deg;
            double Omega = 110.87633970146754/rad2Deg;
            double omega = 285.50181319194064/rad2Deg;
            double f = 0.3759584661778027;

            struct reb_particle flora = reb_tools_orbit_to_particle(r->G, primary, m, a, e, inc, Omega, omega, f);
            reb_add(r,flora);
        }
    }else{
        printf("Found simulation archive. Loaded snapshot at t=%.16f.\n",r->t);
    }

    if (remove(filename) == 0)
      printf("Deleted old sim archive\n");
   else
      printf("Did not find old sim archive\n");

    reb_simulationarchive_automate_interval(r,filename,5e3*sec2year);
    printf("Starting sim at t: %f\n",r->t);
    reb_integrate(r, tmax);
    printf("Final time: %f\n",r->t);
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
