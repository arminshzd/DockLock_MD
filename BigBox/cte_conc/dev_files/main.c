#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "./funcs.h"

// ----------------------------------------------------------------------------------------------------------------------------------------------
// Constants
double temperature = 500;  // K
double kb = 1.38065;  // A^2.g/particle/K/s^2
double Na = 6.022e23;
double step_size = 1e-15 ; // s
double epsilon = 4e-8;  // friction g/s
unsigned long long num_runs = 1e5;
int steps_per_run = 10000;

char filename[] = "output.lammpstrj";
int pbc_on = 1;
double box_size = 60.0; //side/rd;
int n_mols = 60;
int max_growth = 30;

// Umbrella sampling
int US_on = 0;
double krU = 0;
double r0U = 0;
double kqU = 0;
double q0U = 0;
char cfilename[] = "colvars.traj";

// ----------------------------------------------------------------------------------------------------------------------------------------------
// Main Body
int main(){
    pid_t pid = getpid();
    gsl_rng *randgen = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(randgen, pid);
    //gsl_rng_env_setup();
    double m[n_mols*2];  // Mass of particles
    for (int i=0; i < n_mols*2; i++){
        m[i] = 10/Na;
    }
    double beta = 1/temperature/kb;
    double *mp = &m[0];
    double D = 1/beta/epsilon;
    double tau = step_size*D/rd/rd; // dimless parameter
    
    // Init
    int L_o; // Size of the fibril
    int L_n; // Size of the fibril
    double xi[(n_mols+max_growth)*2];
    double yi[(n_mols+max_growth)*2];
    double zi[(n_mols+max_growth)*2];

    //init_pos(xi, yi, zi);
    int inp_status = read_restart("init.coord", xi, yi, zi);
    if (inp_status){
        return 1;
    }
    
    struct Cluster clusters[n_mols];
    struct Cluster *clustersp = &clusters[0];
    int cluster_status = cluster_check(xi, yi, zi, clustersp);
    L_o = fib_size(clustersp, cluster_status);
    // Main Loop
    double potential_energy;
    for (unsigned long long k = 0; k<num_runs; k++){
        for(int i=0;i<steps_per_run;i++){
            // Energy calcs
            potential_energy = u_tot(xi, yi, zi);
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
            // Output screen
            if(i%1000==0){
                FILE *loghandle = fopen("run.log", "a");
                fprintf(loghandle, "%lld\t%f\n", k*steps_per_run+i, potential_energy);
		        fclose(loghandle);
                write_output(filename, k*steps_per_run+i, xi, yi, zi, mp);
	        }
        }
        struct Cluster clusters[n_mols];
        struct Cluster *clustersp = &clusters[0];
        int cluster_status = cluster_check(xi, yi, zi, clustersp);
        L_n = fib_size(clustersp, cluster_status);
        write_colvars_L("colvars.traj", (k+1)*steps_per_run, L_n);
        if (L_n != L_o){
            write_output("finalframe_pre_add.lammpstrj", (k+1)*steps_per_run, xi, yi, zi, mp);
            if (L_n > L_o){
                for (int j=0; j<L_n-L_o; j++){
                    init_pos_mono(xi, yi, zi, n_mols+j);
                }
            }
            n_mols = n_mols + L_n - L_o;
            write_output("finalframe_post_add.lammpstrj", (k+1)*steps_per_run, xi, yi, zi, mp);
            L_o = L_n;
        }
        if (cluster_status != 0){
                write_output("finalframe_pre_break.lammpstrj", (k+1)*steps_per_run, xi, yi, zi, mp);
                break_clusters(xi, yi, zi, clustersp, cluster_status);
                printf("Current fibril size: %d. Resolving %d cluster(s)\n", L_n, cluster_status-1);
                write_output("finalframe_post_break.lammpstrj", (k+1)*steps_per_run, xi, yi, zi, mp);
            }
    }
    return 0;
}

