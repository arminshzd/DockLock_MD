#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "funcs.h"

// ----------------------------------------------------------------------------------------------------------------------------------------------
    // Run Parameters
    double temperature = 500;  // Kelvin
    double kb = 1.38065;  // A^2.g/particle/K/s^2
    double Na = 6.022e23;
    double step_size = 1e-15 ; // s
    double epsilon = 4e-8;  // friction g/s
    char filename[] = "output.lammpstrj";
    char cfilename[] = "colvars.traj";
    double term_dist = 20;
    
    int US_on = 0;
    double krU = 0;
    double r0U = 0;
    double kqU = 0;
    double q0U = 0;

    int pbc_on = 0;
    double box_size = 100; //side/rd;
    int n_mols = 5;

// ----------------------------------------------------------------------------------------------------------------------------------------------
// Main loop
double run(unsigned long long num_steps, int run_num){
    // ----------------------------------------------------------------------------------------------------------------------------------------------
    pid_t pid = getpid();
    gsl_rng *randgen = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(randgen, pid);
    double success = 0;
    double failed = 0;
    double q_fin;
    double r_fin;
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
    double xi[n_mols*2];
    double yi[n_mols*2];
    double zi[n_mols*2];
    
    int inp_status = read_restart("init.coord", xi, yi, zi);
    // The first two are for the monomer and the rest is the fibril
    if (inp_status){
        return 1;
    }

    int r_length = (n_mols-1); // size of the unrolled upper triangle of the distance matrix. look below.
    double r[r_length]; //1-2, 1-3, 1-4, ..., 2-3, 2-4, ..., n-n+1 ..., n-2*n_mols, ...
    double r_min;
    double q;
    // Energy calcs
    for(unsigned long long i=0;i<num_steps;i++){
        double potential_energy = u_tot(xi, yi, zi);
        double kinetic_energy = 0;
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
        //calc_r(xi, yi, zi, r);
        //r_min = arr_min(r, r_length);
        //int max_ind = arr_where(r, r_length, r_min) + 1;
        double comdif[3] = {-((xi[0]+xi[1])-(xi[2]+xi[3]))/2, -((yi[0]+yi[1])-(yi[2]+yi[3]))/2, -((zi[0]+zi[1])-(zi[2]+zi[3]))/2};
        //double ncomdif[3] = {0, 0, 0};
        //mic(comdif, ncomdif);
        r_min = vec_size(comdif);
        double v1[3] = {(xi[0]-xi[1]), (yi[0]-yi[1]), (zi[0]-zi[1])};
        double v2[3] = {(xi[2]-xi[3]), (yi[2]-yi[3]), (zi[2]-zi[3])};
        double cross[3];
        vcross(v1, v2, cross);
        q = vdot(cross, comdif)/r0b/r0b/r_min;
        //char *col_name;
	//sprintf(col_name, "./trajs/traj_%d", run_num);
	//if(i % 25000 == 0){
	//	write_colvars(col_name, i, r_min, q);
	//}
        if(r_min*rd > term_dist || r_min*rd < 1.6){
            printf("%f\t%f\n", q, r_min);
            break;
        }
    }
    //if(r_min*rd < 1.6){
    q_fin = q;
    r_fin = r_min*rd;
    //}
    //else if(r_min*rd > term_dist){
    //        q_fin = 0;
    //        r_fin = 0;
    //}else{
    //    q_fin = -10;
    //    r_fin = -10;
    //}
    FILE *ratehandle = fopen("rate_output.traj", "a");
    fprintf(ratehandle, "%f\t%f\n", q_fin, r_fin);
    fclose(ratehandle);
    return q_fin;
}

int main(){
    return 0;
}
