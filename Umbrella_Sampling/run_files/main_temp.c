#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../../run_files/funcs.h"
// ----------------------------------------------------------------------------------------------------------------------------------------------
// Constants
double temperature = 500;  // Kelvin
double kb = 1.38065;  // A^2.g/particle/K/s^2
double Na = 6.022e23;
double step_size = 1e-15 ; // s
double epsilon = 4e-8;  // friction g/s
unsigned long long num_runs = 1000;
int steps_per_run = 1000;

char filename[] = "output.lammpstrj";
int pbc_on = 0;
double box_size = 100; //side/rd;
int n_mols = 7;

// Umbrella sampling
int US_on = 1;
double krU = 200;
double r0U = 4.0;
double kqU = 100;
double q0U = 1.1;
char cfilename[] = "colvars.traj";
// ----------------------------------------------------------------------------------------------------------------------------------------------
// Main Body

int main(){
    gsl_rng *randgen = gsl_rng_alloc(gsl_rng_taus);
    //gsl_rng_env_setup();
    double m[n_mols*2];  // Mass of particles
    double velx[2*n_mols]; //Velocity
    double vely[2*n_mols];
    double velz[2*n_mols];
    for (int i=0; i < n_mols*2; i++){
        m[i] = 10/Na;
        velx[i] = 0;
        vely[i] = 0;
        velz[i] = 0;
    }

    double beta = 1/temperature/kb;
    double *mp = &m[0];
    double D = 1/beta/epsilon;
    double tau = step_size*D/rd/rd; // dimless parameter
    
    // Init
    int L; // Size of the fibril
    double xi[n_mols*2];
    double yi[n_mols*2];
    double zi[n_mols*2];
    init_pos(xi, yi, zi);
    int inp_status = read_restart("init.coord", xi, yi, zi);
    if (inp_status){
        return 1;
    }
    
    //write_output(filename, 0, xi, yi, zi, m);
    // Main Loop
    double potential_energy;
    double kinetic_energy;
    for (unsigned long long k = 0; k<num_runs; k++){
        for(int i=0;i<steps_per_run;i++){
            // Energy calcs
            potential_energy = u_tot(xi, yi, zi);
            kinetic_energy = 0;
            for (int j=0; j<2*n_mols; j++){
                double vel_vec[3] = {velx[j], vely[j], velz[j]};
                double vel = vec_size(vel_vec);
                kinetic_energy = kinetic_energy + 0.5*m[j]*vel*vel;
            }
            double e_total = potential_energy + kinetic_energy;
            // Backing up pos
            double xo[n_mols*2];
            double yo[n_mols*2];
            double zo[n_mols*2];
            for(int j=0;j<2*n_mols;j++){
                    xo[j] = xi[j];
                    yo[j] = yi[j];
                    zo[j] = zi[j];
            }
            // Integrate one step
            langevin_integrate(xi, yi, zi, tau, randgen);
            // Velocity calculation
            for(int j=0;j<2*n_mols;j++){
                velx[j] = (xi[j] - xo[j])*rd/step_size;
                vely[j] = (yi[j] - yo[j])*rd/step_size;
                velz[j] = (zi[j] - zo[j])*rd/step_size;
            }
            // Output screen
            if(i%100==0){
                //write_output(filename, k*steps_per_run+i, xi, yi, zi, mp);
		if (US_on){
                	double comdif[3] = {-((xi[0]+xi[1])-(xi[2]+xi[3]))/2, -((yi[0]+yi[1])-(yi[2]+yi[3]))/2, -((zi[0]+zi[1])-(zi[2]+zi[3]))/2};
                	double ncomdif[3] = {0, 0, 0};
                	mic(comdif, ncomdif);
                	double r = vec_size(ncomdif);
                	double v1[3] = {(xi[0]-xi[1]), (yi[0]-yi[1]), (zi[0]-zi[1])};
                	double v2[3] = {(xi[2]-xi[3]), (yi[2]-yi[3]), (zi[2]-zi[3])};
                	double nv1[3] = {0, 0, 0};
                	mic(v1, nv1);
                	double nv2[3] = {0, 0, 0};
                	mic(v2, nv2);
                	double cross[3] = {0, 0, 0};
                	vcross(nv1, nv2, cross);
                	double q = vdot(cross, ncomdif)/r0b/r0b/r;
                	//printf("%lld\t%f\t%f\t%f\t%f\n", k*steps_per_run+i, potential_energy, 2*kinetic_energy/3/(n_mols*2)/kb, r, q);
			write_colvars("colvars.traj", (k+1)*steps_per_run+i, r, q);
			}
		}
        }
       // struct Cluster clusters[n_mols];
       // struct Cluster *clustersp = &clusters[0];
       // int cluster_status = cluster_check(xi, yi, zi, clustersp);
       // if (cluster_status != 0){
       //     printf("%d Cluster(s) detected, resolving...\n", cluster_status-1);
       //     write_output("finalframe.lammpstrj", (k+1)*steps_per_run, xi, yi, zi, mp);
       //     L = break_clusters(xi, yi, zi, clustersp, cluster_status);
       //     write_colvars("colvars.traj", (k+1)*steps_per_run, L);
       // }
    }
    return 0;
}

