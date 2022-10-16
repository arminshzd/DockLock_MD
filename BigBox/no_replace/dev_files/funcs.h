#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

struct Cluster {
        int size;
        int members[1000]; // be careful with this. might need to increase it.
};

// ----------------------------------------------------------------------------------------------------------------------------------------------
// Constants
double rd = 2.0;  // length scale of the system A

// LJ potential
double elj = 3.0;
double cofflj = 2.5;
double siglj = 1.0; //2/rd

// COM-LJ potential
double eljc = 5.0; //15
double coffljc = 1.0;
double sigljc = 0.7; //0.7*siglj

// WCA potential
double ewca = 15.0;
double sigwca = 1.4; //1.4*siglj

// Bonded potential
double eb = 1000.0;
double r0b = 1.0; //2/rd

// The Wall
double ea = 6.0;  //from 10.0
double r0a = 1.1;
double sigra = 0.2;
double q0a = -1.0;
double sigqa = 1.2;

// The Channel
double ec = -8.0;
double r0c = 0.7;
double sigrc = 0.5;
double q0c = 0.65;
double sigqc = 0.6;

// Locked state tilting
double etl = -35.0;
double r0t = 0.7;
double sigrt = 0.2;
double q0tl = -1.0;
double sigqt = 0.6;


extern int pbc_on;
extern double box_size;
extern int n_mols;

// Umbrella sampling
extern int US_on;
extern double krU;
extern double r0U;
extern double kqU;
extern double q0U;

// ----------------------------------------------------------------------------------------------------------------------------------------------
// Utility functions

double vec_size(double *vec){
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

void vcross(double *v1, double *v2, double *res){
    res[0] = v1[1]*v2[2] - v1[2]*v2[1];
    res[1] = v1[2]*v2[0] - v1[0]*v2[2];
    res[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

double vdot(double *v1, double *v2){
    double res = 0;
    for(int i=0; i<3; i++) {
        res += v1[i]*v2[i];
    }
    return res;
}

// double f_norm(double r, double mu, double sig){
// 	return (1/sqrt(2*M_PI*sig*sig))*exp(-(r-mu)*(r-mu)/(2*sig*sig));
// }

int mic(double *vec, double *resvec){
    if (pbc_on){
    double x_rsize = 1.0/box_size;
    for (int i=0; i<3; i++){
        resvec[i] = vec[i] - box_size*round(vec[i]*x_rsize);
    }
    } else {
        for (int i=0; i<3; i++){
        resvec[i] = vec[i];
        }
    }
    return 0;
}

double pbc_cord(double x){
    double xn;
    if (x > box_size/2.0){
        xn = x - box_size;
    }else if(x <= -box_size/2.0){
        xn = x + box_size;
    }else{
        xn = x;
    }
    return xn;
}

int valueinarray(int val, int *arr, int arrsize){
    //for(int i = 0; i < sizeof(arr) / sizeof(arr[0]); i++){
    for(int i=0; i< arrsize; i++){
        if(arr[i] == val){
            return 1;
        }
    }
    return 0;
}

int init_pos(double *x, double *y, double *z){
    gsl_rng *randgen = gsl_rng_alloc(gsl_rng_taus);
    double xc[n_mols];
    double yc[n_mols];
    double zc[n_mols];
    double dirx[n_mols];
    double diry[n_mols];
    double dirz[n_mols];
    double cnt = 0;
    int len = 0;
    int flag1 = 1;
    int flag2 = 1;
    while (flag1){
        double rndx = box_size*gsl_rng_uniform(randgen)*pow(-1, round(gsl_rng_uniform(randgen)))/2.0;
        double rndy = box_size*gsl_rng_uniform(randgen)*pow(-1, round(gsl_rng_uniform(randgen)))/2.0;
        double rndz = box_size*gsl_rng_uniform(randgen)*pow(-1, round(gsl_rng_uniform(randgen)))/2.0;
        if (len != 0){
            for (int i = 0; i < len; i++){
                double vec1[3] = {xc[i] - rndx, yc[i] - rndy, zc[i] - rndz};
                if (vec_size(vec1) <= 1){
                    flag2 = 0;
                    break;
                }
            }
            if (flag2){
                xc[len] = rndx;
                yc[len] = rndy;
                zc[len] = rndz;
                len = len + 1;
            }
            if (len == n_mols){
                flag1 = 0;
            }
        
        } else {
            len = len + 1;
            xc[0] = rndx;
            yc[0] = rndy;
            zc[0] = rndz;
        }
    }
    for (int i=0; i < n_mols; i++){
        double sinth = pow(-1, round(gsl_rng_uniform(randgen)))*gsl_rng_uniform(randgen);
        double costh = sqrt(1-pow(sinth,2))*pow(-1, round(gsl_rng_uniform(randgen)));
        double sinphi = pow(-1, round(gsl_rng_uniform(randgen)))*gsl_rng_uniform(randgen);
        double cosphi = sqrt(1-pow(sinphi,2))*pow(-1, round(gsl_rng_uniform(randgen)));
        double rvecx = sinphi*costh/2.0;
        double rvecy = sinphi*sinth/2.0;
        double rvecz = cosphi/2.0;
        x[2*i] = xc[i] + rvecx;
        x[2*i+1] = xc[i] - rvecx;
        y[2*i] = yc[i] + rvecy;
        y[2*i+1] = yc[i] - rvecy;
        z[2*i] = zc[i] + rvecz;
        z[2*i+1] = zc[i] - rvecz;
    }
    return 0;
}

void init_pos_mono(double *x, double *y, double *z, int mono_id){
    gsl_rng *randgen = gsl_rng_alloc(gsl_rng_taus);
    double xc;
    double yc;
    double zc;
    double dirx;
    double diry;
    double dirz;
    double cnt = 0;
    int len = 0;
    int flag = 1;
    while (flag){
        double rndx = box_size*gsl_rng_uniform(randgen)*pow(-1, round(gsl_rng_uniform(randgen)))/2.0;
        double rndy = box_size*gsl_rng_uniform(randgen)*pow(-1, round(gsl_rng_uniform(randgen)))/2.0;
        double rndz = box_size*gsl_rng_uniform(randgen)*pow(-1, round(gsl_rng_uniform(randgen)))/2.0;
        flag = 0;
        for (int i = 0; i < n_mols*2; i++){
            double difvec[3] = {x[i] - rndx, y[i] - rndy, z[i] - rndz};
            double difvecm[3] = {0, 0, 0};
            mic(difvec, difvecm);
            if (vec_size(difvecm) <= 5){
                flag = 1;
                break;
            }
        }
        xc = rndx;
        yc = rndy;
        zc = rndz;
    }
    double sinth = pow(-1, round(gsl_rng_uniform(randgen)))*gsl_rng_uniform(randgen);
    double costh = sqrt(1-pow(sinth,2))*pow(-1, round(gsl_rng_uniform(randgen)));
    double sinphi = pow(-1, round(gsl_rng_uniform(randgen)))*gsl_rng_uniform(randgen);
    double cosphi = sqrt(1-pow(sinphi,2))*pow(-1, round(gsl_rng_uniform(randgen)));
    double rvecx = sinphi*costh/2.0;
    double rvecy = sinphi*sinth/2.0;
    double rvecz = cosphi/2.0;
    double vectest[3] = {rvecx, rvecy, rvecz};
    x[2*mono_id] = xc + rvecx;
    x[2*mono_id+1] = xc - rvecx;
    y[2*mono_id] = yc + rvecy;
    y[2*mono_id+1] = yc - rvecy;
    z[2*mono_id] = zc + rvecz;
    z[2*mono_id+1] = zc - rvecz;
}


void calc_rq(double *xi, double *yi, double *zi, double *r, double *q){
    int r_cnt = 0;
    for (int i=0; i < 1; i++){ // EDITED TO GIVE DISTANCES FROM THE FIRST MOLECULE ONLY
        for (int j=i+1; j < n_mols; j++){
            if (i != j){
                double tempx[4] = {xi[2*i], xi[2*i+1], xi[2*j], xi[2*j+1]};
                double tempy[4] = {yi[2*i], yi[2*i+1], yi[2*j], yi[2*j+1]};
                double tempz[4] = {zi[2*i], zi[2*i+1], zi[2*j], zi[2*j+1]};
                double comdif[3] = {-((tempx[0]+tempx[1])-(tempx[2]+tempx[3]))/2.0, -((tempy[0]+tempy[1])-(tempy[2]+tempy[3]))/2.0, -((tempz[0]+tempz[1])-(tempz[2]+tempz[3]))/2.0};
                double ncomdif[3] = {0, 0, 0};
                mic(comdif, ncomdif);
                r[r_cnt] = vec_size(ncomdif);
                double v1[3] = {(tempx[0]-tempx[1]), (tempy[0]-tempy[1]), (tempz[0]-tempz[1])};
                double v2[3] = {(tempx[2]-tempx[3]), (tempy[2]-tempy[3]), (tempz[2]-tempz[3])};
                double nv1[3] = {0, 0, 0};
                mic(v1, nv1);
                double nv2[3] = {0, 0, 0};
                mic(v2, nv2);
                double cross[3] = {0, 0, 0};
                vcross(nv1, nv2, cross);
                q[r_cnt] = vdot(cross, ncomdif)/r0b/r0b/r[r_cnt];
                r_cnt++;
            }
        }
    }
}

double arr_min(double *arr, int arrlen){
    double min = arr[0];
    for (int i=1; i<arrlen; i++){
        if (min > arr[i]){
            min = arr[i];
        }
    }
    return min;
}

int arr_where(double *arr, int arrlen, double val){
    for (int i=0; i<arrlen; i++){
        if (arr[i] == val){
            return i;
        }
    }
}

int cluster_check(double *x, double *y, double *z, struct Cluster *clustersp){
    void calc_rq(double *xi, double *yi, double *zi, double *r, double *q){
        int r_cnt = 0;
        for (int i=0; i < n_mols; i++){
            for (int j=i+1; j < n_mols; j++){
                if (i != j){
                    double tempx[4] = {xi[2*i], xi[2*i+1], xi[2*j], xi[2*j+1]};
                    double tempy[4] = {yi[2*i], yi[2*i+1], yi[2*j], yi[2*j+1]};
                    double tempz[4] = {zi[2*i], zi[2*i+1], zi[2*j], zi[2*j+1]};
                    double comdif[3] = {-((tempx[0]+tempx[1])-(tempx[2]+tempx[3]))/2.0, -((tempy[0]+tempy[1])-(tempy[2]+tempy[3]))/2.0, -((tempz[0]+tempz[1])-(tempz[2]+tempz[3]))/2.0};
                    double ncomdif[3] = {0, 0, 0};
                    mic(comdif, ncomdif);
                    r[r_cnt] = vec_size(ncomdif);
                    double v1[3] = {(tempx[0]-tempx[1]), (tempy[0]-tempy[1]), (tempz[0]-tempz[1])};
                    double v2[3] = {(tempx[2]-tempx[3]), (tempy[2]-tempy[3]), (tempz[2]-tempz[3])};
                    double nv1[3] = {0, 0, 0};
                    mic(v1, nv1);
                    double nv2[3] = {0, 0, 0};
                    mic(v2, nv2);
                    double cross[3] = {0, 0, 0};
                    vcross(nv1, nv2, cross);
                    q[r_cnt] = vdot(cross, ncomdif)/r0b/r0b/r[r_cnt];
                    r_cnt++;
                }
            }
        }
    }
    void neighbor_list(double *r, double *q, int *cntlist, int rlen){
        int ii = 0;
        int jj = ii+1;
        for (int i=0; i<rlen; i++){
            if (r[i] < 1.50 && fabs(q[i])>0.70){
                cntlist[ii*n_mols+jj] = 1;
            }
            if (jj<n_mols-1){
                jj += 1;
            } else {
                ii += 1;
                jj = ii + 1;
            }
        }
    }

    void find_clst(int *nlist, int ind, struct Cluster *clstr){
        int row_tracker = ind; // tracking position in the unrolled neighbor list
        clstr->members[0] = ind; // adding the initial monomer to the cluster list
        clstr->size = 1; // initiating the cluster size
        int aflag = 1; // member addition flag
        int cflag = 1; // member checking flag
        int ccount = 0; // position in members list
        while (aflag || cflag){
            aflag = 0;
            cflag = 0;
            for (int i=0; i<n_mols; i++){
                //printf("%d\n", nlist[i + row_tracker*n_mols]);
                if (nlist[i + row_tracker*n_mols]){
                    if (!valueinarray(i, clstr->members, clstr->size)){
                        aflag = 1;
                        clstr->members[clstr->size] = i;
                        clstr->size += 1;
                        //printf("%d\n", clstr->size);
                    }
                }
            }
            if (ccount < (clstr->size)-1){
                cflag = 1;
                row_tracker = clstr->members[ccount+1];
                ccount += 1;
            }
            // for (int i = 0; i<clstr->size; i++){
            //     printf("%d\t", clstr->members[i]);
            // }
            // printf("\n");
            // printf("%d\n", ccount);
            // printf("%d\n", clstr->size);
            // printf("%d\n", valueinarray(4, clstr->members, clstr->size));
        }
    }

    int cluster_list(int *nlist, struct Cluster *clstrs){
        int ind = 0;
        int clstrno = 0;
        int asgn_list[n_mols]; // assignment list. assigns each molecule to a cluster (1, 2, 3, ...) or -1 if its a monomer.
        for (int i=0; i<n_mols; i++){
            asgn_list[i] = -1;
        }
        while (ind < n_mols-1){
            int sum = 0;
            for (int i=0; i<n_mols; i++){
                sum += nlist[ind*n_mols+i];
            }
            while(ind < n_mols-1 && sum == 0){
                ind += 1;
                sum = 0;
                for (int i=0; i<n_mols; i++){
                    sum += nlist[ind*n_mols+i];
                }
            }
            if (asgn_list[ind] == -1 && sum != 0){
                struct Cluster cclstr;
                struct Cluster *cclstrp = &cclstr;
                find_clst(nlist, ind, cclstrp);
                clstrs[clstrno] = cclstr;
                clstrno += 1;
                for (int i=0; i<cclstr.size; i++){
                    asgn_list[cclstr.members[i]] = clstrno;
                }
                
            }
            ind += 1;
        }
        return clstrno;
    }
    
    int r_length = (n_mols-1)*n_mols/2; // size of the unrolled upper triangle of the distance matrix. look below.
    double r[r_length]; //1-2, 1-3, 1-4, ..., 2-3, 2-4, ..., n-n+1 ..., n-2*n_mols, ...
    double q[r_length]; //1-2, 1-3, 1-4, ..., 2-3, 2-4, ..., n-n+1 ..., n-2*n_mols, ...
    calc_rq(x, y, z, r, q);
    int cntlist[n_mols][n_mols]; // neighbors matrix, 1 if they are, 0 if they're not
    for (int i=0; i<n_mols; i++){
        for (int j=0; j<n_mols; j++){
            cntlist[i][j] = 0; // setting intial values to 0
        }
    }
    int *cntlistp = &cntlist[0][0];
    neighbor_list(r, q, cntlistp, r_length); // finding the neighbors. populates the upper triangle only
    for (int i=0; i<n_mols; i++){
        for (int j=0; j<n_mols; j++){
            if (i != j){
                cntlist[j][i] = cntlist[i][j]; // making the matrix full (from the upper triangle part)
            }
        }
    }

    // for (int i=0; i<n_mols; i++){
    //     for (int j=0; j<n_mols; j++){
    //         printf("%d\t", cntlist[i][j]);
    //     }
    //     printf("\n");
    // }

    int no_clusters = cluster_list(cntlistp, clustersp); // PROBLEM HERE
    if (no_clusters > 1){
        return no_clusters;
    }
    return 0;
}

void fix_cluster(double *x, double *y, double *z, struct Cluster *clusterp){
    for (int i=0; i<clusterp->size; i++){
        init_pos_mono(x, y, z, clusterp->members[i]);
    }
}

int break_clusters(double *x, double *y, double *z, struct Cluster *clusters, int n_clusters){
    int max_size = clusters[0].size;
    for (int i=0; i<n_clusters; i++){ 
        if (max_size < clusters[i].size){
            max_size = clusters[i].size;
        }
    }
    int multi_max = 0;
    for (int i=0; i<n_clusters; i++){ 
        if (max_size == clusters[i].size){
            multi_max += 1;
        }
    }
    for (int i=0; i<n_clusters; i++){
        if (clusters[i].size != max_size){
            struct Cluster *clstrp = &clusters[i];
            fix_cluster(x, y, z, clstrp);
        } else if (multi_max > 1) {
            multi_max -= 1;
            struct Cluster *clstrp = &clusters[i];
            fix_cluster(x, y, z, clstrp);
        }
    }
    return max_size;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------
// Potentials
    //Umbrella Sampling (ONLY FOR 2 MOLECULES!!!!!)
double U_U_distance(double *x, double *y, double *z, double k, double r0){
    double comdif[3] = {-((x[0]+x[1])-(x[2]+x[3]))/2.0, -((y[0]+y[1])-(y[2]+y[3]))/2.0, -((z[0]+z[1])-(z[2]+z[3]))/2.0};
    double ncomdif[3] = {0, 0, 0};
    mic(comdif, ncomdif);
    return 0.5*k*(vec_size(ncomdif)-r0)*(vec_size(ncomdif)-r0);
}
    //Umbrella Sampling (ONLY FOR 2 MOLECULES!!!!!)
double U_U_q(double *x, double *y, double *z, double k, double q0){
    double comdif[3] = {-((x[0]+x[1])-(x[2]+x[3]))/2.0, -((y[0]+y[1])-(y[2]+y[3]))/2.0, -((z[0]+z[1])-(z[2]+z[3]))/2.0};
    double ncomdif[3] = {0, 0, 0};
    mic(comdif, ncomdif);
    double r = vec_size(ncomdif);
    double v1[3] = {(x[0]-x[1]), (y[0]-y[1]), (z[0]-z[1])};
    double v2[3] = {(x[2]-x[3]), (y[2]-y[3]), (z[2]-z[3])};
    double nv1[3] = {0, 0, 0};
    mic(v1, nv1);
    double nv2[3] = {0, 0, 0};
    mic(v2, nv2);
    double cross[3] = {0, 0, 0};
    vcross(nv1, nv2, cross);
    double q = vdot(cross, ncomdif)/r0b/r0b/r;
    return 0.5*k*(q-q0)*(q-q0);
}

double uljmain(double rd, double e, double sig){
    return 4*e*(pow(((sig/pow(2,1.0/6.0))/rd),12) - pow(((sig/pow(2,1.0/6.0))/rd),6));
}

double ulj(double *v1, double *v2){
    double distv[3] = {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
    double ndistv[3] = {0, 0, 0};
    mic(distv, ndistv);
    double dist = vec_size(ndistv);
    if (dist <= cofflj){
        return uljmain(dist, elj, siglj) - uljmain(cofflj, elj, siglj);
    } else{
        return 0;
    }
}

double uwca(double *v1, double *v2){
    double distv[3] = {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
    double ndistv[3] = {0, 0, 0};
    mic(distv, ndistv);
    double dist = vec_size(ndistv);
    if (dist <= sigwca){
        return uljmain(dist, ewca, sigwca) - uljmain(sigwca, ewca, sigwca);
    } else{
        return 0;
    }
}

double uljc(double *x, double *y, double *z){
    double comdif[3] = {-((x[0]+x[1])-(x[2]+x[3]))/2.0, -((y[0]+y[1])-(y[2]+y[3]))/2.0, -((z[0]+z[1])-(z[2]+z[3]))/2.0};
    double ncomdif[3] = {0, 0, 0};
    mic(comdif, ncomdif);
    double r = vec_size(ncomdif);
    if (r <= coffljc){
        return uljmain(r, eljc, sigljc) - uljmain(coffljc, eljc, sigljc);
    } else{
        return 0;
    }
}

double ubonded(double *v1, double *v2){
    double distv[3] = {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
    double ndistv[3] = {0, 0, 0};
    mic(distv, ndistv);
    double dist = vec_size(ndistv);
    return 0.5*eb*pow((dist-r0b),2);
}

double uangle(double *x, double *y, double *z, double e, double r0, double sr, double q0, double sq){
    double comdif[3] = {-((x[0]+x[1])-(x[2]+x[3]))/2.0, -((y[0]+y[1])-(y[2]+y[3]))/2.0, -((z[0]+z[1])-(z[2]+z[3]))/2.0};
    double ncomdif[3] = {0, 0, 0};
    mic(comdif, ncomdif);
    double r = vec_size(ncomdif);
    double v1[3] = {(x[0]-x[1]), (y[0]-y[1]), (z[0]-z[1])};
    double v2[3] = {(x[2]-x[3]), (y[2]-y[3]), (z[2]-z[3])};
    double nv1[3] = {0, 0, 0};
    mic(v1, nv1);
    double nv2[3] = {0, 0, 0};
    mic(v2, nv2);
    double cross[3] = {0, 0, 0};
    vcross(nv1, nv2, cross);
    double q = vdot(cross, ncomdif)/r0b/r0b/r;
    double f_com = exp(-(r-r0)*(r-r0)/(2*sr*sr));
    double f_q = exp(-(q-q0)*(q-q0)/(2*sq*sq));
    return e*f_com*f_q;
}

double u_tot(double *x, double *y, double *z){
    double bonded = 0;
    for (int i=0; i < 2*n_mols-1; i+=2){
        double vb1[3] = {x[i], y[i], z[i]};
        double vb2[3] = {x[i+1], y[i+1], z[i+1]};
        bonded += ubonded(vb1, vb2);
    }

    double non_bonded = 0;
    for (int i=0; i < 2*n_mols; i++){
        for (int j=0; j < 2*n_mols; j++){
            if (i != j){
                if (i%4==0){
                    if (j%4==2 || j%4==3){
                        double vlj1[3] = {x[i], y[i], z[i]};
                        double vlj2[3] = {x[j], y[j], z[j]};
                        non_bonded = non_bonded + ulj(vlj1, vlj2);
                    } else if(j!=i+1){
                        double vlj1[3] = {x[i], y[i], z[i]};
                        double vlj2[3] = {x[j], y[j], z[j]};
                        non_bonded = non_bonded + uwca(vlj1, vlj2);
                    }
                }else if(i%4==1){
                    if (j%4==2 || j%4==3){
                        double vlj1[3] = {x[i], y[i], z[i]};
                        double vlj2[3] = {x[j], y[j], z[j]};
                        non_bonded = non_bonded + ulj(vlj1, vlj2);
                    } else if(j!=i-1){
                        double vlj1[3] = {x[i], y[i], z[i]};
                        double vlj2[3] = {x[j], y[j], z[j]};
                        non_bonded = non_bonded + uwca(vlj1, vlj2);
                    }
                }else if(i%4==2){
                    if (j%4==0 || j%4==1){
                        double vlj1[3] = {x[i], y[i], z[i]};
                        double vlj2[3] = {x[j], y[j], z[j]};
                        non_bonded = non_bonded + ulj(vlj1, vlj2);
                    } else if(j!=i+1){
                        double vlj1[3] = {x[i], y[i], z[i]};
                        double vlj2[3] = {x[j], y[j], z[j]};
                        non_bonded = non_bonded + uwca(vlj1, vlj2);
                    }
                }else if(i%4==3){
                    if (j%4==0 || j%4==1){
                        double vlj1[3] = {x[i], y[i], z[i]};
                        double vlj2[3] = {x[j], y[j], z[j]};
                        non_bonded = non_bonded + ulj(vlj1, vlj2);
                    } else if(j!=i-1){
                        double vlj1[3] = {x[i], y[i], z[i]};
                        double vlj2[3] = {x[j], y[j], z[j]};
                        non_bonded = non_bonded + uwca(vlj1, vlj2);
                    }
                }else{
                    printf("Something went wrong...");
                }
            }
        }
    }
    double angle = 0;
    for (int i=0; i < n_mols; i++){
        for (int j=0; j < n_mols; j++){
            if (i != j && (2*i)%4 != (2*j)%4){
                double tempx[4] = {x[2*i], x[2*i+1], x[2*j], x[2*j+1]};
                double tempy[4] = {y[2*i], y[2*i+1], y[2*j], y[2*j+1]};
                double tempz[4] = {z[2*i], z[2*i+1], z[2*j], z[2*j+1]};
                angle = angle + uangle(tempx, tempy, tempz, ea, r0a, sigra, q0a, sigqa); // The Wall
		angle = angle + uljc(tempx, tempy, tempz); // C-to-C LJ
                angle = angle + uangle(tempx, tempy, tempz, ec, r0c, sigrc, q0c, sigqc); // The Channel
                angle = angle + uangle(tempx, tempy, tempz, etl, r0t, sigrt, q0tl, sigqt); // locked state tilting
            }
        }
    }
    double umbrella; // FOR 2 MOLECULES ONLY
    if (US_on){
        umbrella = U_U_distance(x, y, z, krU, r0U) + U_U_q(x, y, z, kqU, q0U);
    }else{
        umbrella = 0;
    }
    // Umbrella sampling
    return umbrella + non_bonded/2.0 + bonded + angle/2.0;
}
// ----------------------------------------------------------------------------------------------------------------------------------------------
// Forces
// Umbrella sampling (2 MOLECULES ONLY)
void f_U_distance(double *x, double *y, double *z, double *rx, double *ry, double *rz, double k, double r0){
    double comdif[3] = {-((x[0]+x[1])-(x[2]+x[3]))/2.0, -((y[0]+y[1])-(y[2]+y[3]))/2.0, -((z[0]+z[1])-(z[2]+z[3]))/2.0};
    double ncomdif[3] = {0, 0, 0};
    mic(comdif, ncomdif);
    double r = vec_size(ncomdif);
    rx[0] = -k * (r-r0) * -1 * ncomdif[0] * (1/r) * 0.5;
    ry[0] = -k * (r-r0) * -1 * ncomdif[1] * (1/r) * 0.5;
    rz[0] = -k * (r-r0) * -1 * ncomdif[2] * (1/r) * 0.5;
    rx[1] = -k * (r-r0) * -1 * ncomdif[0] * (1/r) * 0.5;
    ry[1] = -k * (r-r0) * -1 * ncomdif[1] * (1/r) * 0.5;
    rz[1] = -k * (r-r0) * -1 * ncomdif[2] * (1/r) * 0.5;
    rx[2] = -k * (r-r0) * ncomdif[0] * (1/r) * 0.5;
    ry[2] = -k * (r-r0) * ncomdif[1] * (1/r) * 0.5;
    rz[2] = -k * (r-r0) * ncomdif[2] * (1/r) * 0.5;
    rx[3] = -k * (r-r0) * ncomdif[0] * (1/r) * 0.5;
    ry[3] = -k * (r-r0) * ncomdif[1] * (1/r) * 0.5;
    rz[3] = -k * (r-r0) * ncomdif[2] * (1/r) * 0.5;
}

// Umbrella sampling (2 MOLECULES ONLY)
void f_U_q(double *x, double *y, double *z, double *rx, double *ry, double *rz, double k, double q0){
    double comdif[3] = {-((x[0]+x[1])-(x[2]+x[3]))/2.0, -((y[0]+y[1])-(y[2]+y[3]))/2.0, -((z[0]+z[1])-(z[2]+z[3]))/2.0};
    double ncomdif[3] = {0, 0, 0};
    mic(comdif, ncomdif);
    double r = vec_size(ncomdif);
    double v1[3] = {(x[0]-x[1]), (y[0]-y[1]), (z[0]-z[1])};
    double v2[3] = {(x[2]-x[3]), (y[2]-y[3]), (z[2]-z[3])};
    double nv1[3] = {0, 0, 0};
    mic(v1, nv1);
    double nv2[3] = {0, 0, 0};
    mic(v2, nv2);
    double cross[3] = {0, 0, 0};
    vcross(nv1, nv2, cross);
    double q = vdot(cross, ncomdif)/r0b/r0b/r;
    double dxxdx[4][3] = {{0, nv2[2], -nv2[1]}, {0, -nv2[2], nv2[1]}, {0, -nv1[2], nv1[1]}, {0, nv1[2], -nv1[1]}};
    double dyxdx[4][3] = {{-nv2[2], 0, nv2[0]}, {nv2[2], 0, -nv2[0]}, {nv1[2], 0, -nv1[0]}, {-nv1[2], 0, nv1[0]}};
    double dzxdx[4][3] = {{nv2[1], -nv2[0], 0}, {-nv2[1], nv2[0], 0}, {-nv1[1], nv1[0], 0}, {nv1[1], -nv1[0], 0}};
    double dqdx[4] = {ncomdif[1]*dyxdx[0][0]+ncomdif[2]*dzxdx[0][0]+cross[0]*(-0.5), ncomdif[1]*dyxdx[1][0]+ncomdif[2]*dzxdx[1][0]+cross[0]*(-0.5), ncomdif[1]*dyxdx[2][0]+ncomdif[2]*dzxdx[2][0]+cross[0]*(0.5), ncomdif[1]*dyxdx[3][0]+ncomdif[2]*dzxdx[3][0]+cross[0]*(0.5)};
    double dqdy[4] = {ncomdif[0]*dxxdx[0][1]+ncomdif[2]*dzxdx[0][1]+cross[1]*(-0.5), ncomdif[0]*dxxdx[1][1]+ncomdif[2]*dzxdx[1][1]+cross[1]*(-0.5), ncomdif[0]*dxxdx[2][1]+ncomdif[2]*dzxdx[2][1]+cross[1]*(0.5), ncomdif[0]*dxxdx[3][1]+ncomdif[2]*dzxdx[3][1]+cross[1]*(0.5)};
    double dqdz[4] = {ncomdif[0]*dyxdx[0][2]+ncomdif[1]*dzxdx[0][0]+cross[2]*(-0.5), ncomdif[0]*dyxdx[1][2]+ncomdif[1]*dzxdx[1][0]+cross[2]*(-0.5), ncomdif[0]*dyxdx[2][2]+ncomdif[1]*dzxdx[2][0]+cross[2]*(0.5), ncomdif[0]*dyxdx[3][2]+ncomdif[1]*dzxdx[3][0]+cross[2]*(0.5)};
    rx[0] = -k * (q-q0)*(dqdx[0]/r/r0b/r0b + q/r/r*ncomdif[0]*(1)*0.5);
    ry[0] = -k * (q-q0)*(dqdy[0]/r/r0b/r0b + q/r/r*ncomdif[1]*(1)*0.5);
    rz[0] = -k * (q-q0)*(dqdz[0]/r/r0b/r0b + q/r/r*ncomdif[2]*(1)*0.5);
    rx[1] = -k * (q-q0)*(dqdx[1]/r/r0b/r0b + q/r/r*ncomdif[0]*(1)*0.5);
    ry[1] = -k * (q-q0)*(dqdy[1]/r/r0b/r0b + q/r/r*ncomdif[1]*(1)*0.5);
    rz[1] = -k * (q-q0)*(dqdz[1]/r/r0b/r0b + q/r/r*ncomdif[2]*(1)*0.5);
    rx[2] = -k * (q-q0)*(dqdx[2]/r/r0b/r0b + q/r/r*ncomdif[0]*(-1)*0.5);
    ry[2] = -k * (q-q0)*(dqdy[2]/r/r0b/r0b + q/r/r*ncomdif[1]*(-1)*0.5);
    rz[2] = -k * (q-q0)*(dqdz[2]/r/r0b/r0b + q/r/r*ncomdif[2]*(-1)*0.5);
    rx[3] = -k * (q-q0)*(dqdx[3]/r/r0b/r0b + q/r/r*ncomdif[0]*(-1)*0.5);
    ry[3] = -k * (q-q0)*(dqdy[3]/r/r0b/r0b + q/r/r*ncomdif[1]*(-1)*0.5);
    rz[3] = -k * (q-q0)*(dqdz[3]/r/r0b/r0b + q/r/r*ncomdif[2]*(-1)*0.5);
}

void fangle(double *x, double *y, double *z, double *rx, double *ry, double *rz, double e, double r0, double sr, double q0, double sq){
    double comdif[3] = {-((x[0]+x[1])-(x[2]+x[3]))/2.0, -((y[0]+y[1])-(y[2]+y[3]))/2.0, -((z[0]+z[1])-(z[2]+z[3]))/2.0};
    double ncomdif[3] = {0, 0, 0};
    mic(comdif, ncomdif);
    double r = vec_size(ncomdif);
    double v1[3] = {(x[0]-x[1]), (y[0]-y[1]), (z[0]-z[1])};
    double v2[3] = {(x[2]-x[3]), (y[2]-y[3]), (z[2]-z[3])};
    double nv1[3] = {0, 0, 0};
    mic(v1, nv1);
    double nv2[3] = {0, 0, 0};
    mic(v2, nv2);
    double cross[3] = {0, 0, 0};
    vcross(nv1, nv2, cross);
    double q = vdot(cross, ncomdif)/r0b/r0b/r;
    double f_com = exp(-(r-r0)*(r-r0)/(2*sr*sr));
    double f_com_prime_r = -1.00*f_com*(r-r0)/(sr*sr);
    double f_com_prime_q = 0;
    double f_q = exp(-(q-q0)*(q-q0)/(2*sq*sq));
    double f_q_prime_q = -1.00*f_q*(q-q0)/(sq*sq)/r;
    double f_q_prime_r = f_q*(q-q0)/(sq*sq)*q/r;
    double dUdrc = e*(f_com_prime_r*f_q + f_com*f_q_prime_r);
    double dUdq = e*(f_com_prime_q*f_q + f_com*f_q_prime_q); 
    double dxxdx[4][3] = {{0, nv2[2], -nv2[1]}, {0, -nv2[2], nv2[1]}, {0, -nv1[2], nv1[1]}, {0, nv1[2], -nv1[1]}};
    double dyxdx[4][3] = {{-nv2[2], 0, nv2[0]}, {nv2[2], 0, -nv2[0]}, {nv1[2], 0, -nv1[0]}, {-nv1[2], 0, nv1[0]}};
    double dzxdx[4][3] = {{nv2[1], -nv2[0], 0}, {-nv2[1], nv2[0], 0}, {-nv1[1], nv1[0], 0}, {nv1[1], -nv1[0], 0}};
    double dqdx[4] = {ncomdif[1]*dyxdx[0][0]+ncomdif[2]*dzxdx[0][0]+cross[0]*(-0.5), ncomdif[1]*dyxdx[1][0]+ncomdif[2]*dzxdx[1][0]+cross[0]*(-0.5), ncomdif[1]*dyxdx[2][0]+ncomdif[2]*dzxdx[2][0]+cross[0]*(0.5), ncomdif[1]*dyxdx[3][0]+ncomdif[2]*dzxdx[3][0]+cross[0]*(0.5)};
    double dqdy[4] = {ncomdif[0]*dxxdx[0][1]+ncomdif[2]*dzxdx[0][1]+cross[1]*(-0.5), ncomdif[0]*dxxdx[1][1]+ncomdif[2]*dzxdx[1][1]+cross[1]*(-0.5), ncomdif[0]*dxxdx[2][1]+ncomdif[2]*dzxdx[2][1]+cross[1]*(0.5), ncomdif[0]*dxxdx[3][1]+ncomdif[2]*dzxdx[3][1]+cross[1]*(0.5)};
    double dqdz[4] = {ncomdif[0]*dyxdx[0][2]+ncomdif[1]*dzxdx[0][0]+cross[2]*(-0.5), ncomdif[0]*dyxdx[1][2]+ncomdif[1]*dzxdx[1][0]+cross[2]*(-0.5), ncomdif[0]*dyxdx[2][2]+ncomdif[1]*dzxdx[2][0]+cross[2]*(0.5), ncomdif[0]*dyxdx[3][2]+ncomdif[1]*dzxdx[3][0]+cross[2]*(0.5)};
    rx[0] = -(dUdrc*ncomdif[0]*(-1)*0.5/r + dUdq*dqdx[0]);
    ry[0] = -(dUdrc*ncomdif[1]*(-1)*0.5/r + dUdq*dqdy[0]);
    rz[0] = -(dUdrc*ncomdif[2]*(-1)*0.5/r + dUdq*dqdz[0]);
    rx[1] = -(dUdrc*ncomdif[0]*(-1)*0.5/r + dUdq*dqdx[1]);
    ry[1] = -(dUdrc*ncomdif[1]*(-1)*0.5/r + dUdq*dqdy[1]);
    rz[1] = -(dUdrc*ncomdif[2]*(-1)*0.5/r + dUdq*dqdz[1]);
    rx[2] = -(dUdrc*ncomdif[0]*(1)*0.5/r + dUdq*dqdx[2]);
    ry[2] = -(dUdrc*ncomdif[1]*(1)*0.5/r + dUdq*dqdy[2]);
    rz[2] = -(dUdrc*ncomdif[2]*(1)*0.5/r + dUdq*dqdz[2]);
    rx[3] = -(dUdrc*ncomdif[0]*(1)*0.5/r + dUdq*dqdx[3]);
    ry[3] = -(dUdrc*ncomdif[1]*(1)*0.5/r + dUdq*dqdy[3]);
    rz[3] = -(dUdrc*ncomdif[2]*(1)*0.5/r + dUdq*dqdz[3]);
}

void flj(double *v1, double *v2, double *res){
    double distv[3] = {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
    double ndistv[3] = {0, 0, 0};
    mic(distv, ndistv);
    double dist = vec_size(ndistv);
    if (dist <= cofflj){
        res[0] = -4*elj*ndistv[0]*(6*pow((siglj/pow(2,1.0/6.0)),6)/pow(dist,7) - 12*pow((siglj/pow(2,1.0/6.0)),12)/pow(dist,13))/dist;
        res[1] = -4*elj*ndistv[1]*(6*pow((siglj/pow(2,1.0/6.0)),6)/pow(dist,7) - 12*pow((siglj/pow(2,1.0/6.0)),12)/pow(dist,13))/dist;
        res[2] = -4*elj*ndistv[2]*(6*pow((siglj/pow(2,1.0/6.0)),6)/pow(dist,7) - 12*pow((siglj/pow(2,1.0/6.0)),12)/pow(dist,13))/dist;
    } else{
        res[0] = 0;
        res[1] = 0;
        res[2] = 0;
    }
}

void fwca(double *v1, double *v2, double *res){
    double distv[3] = {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
    double ndistv[3] = {0, 0, 0};
    mic(distv, ndistv);
    double dist = vec_size(ndistv);
    if (dist <= sigwca){
        res[0] = -4*ewca*ndistv[0]*(6*pow((sigwca/pow(2,1.0/6.0)),6)/pow(dist,7) - 12*pow((sigwca/pow(2,1.0/6.0)),12)/pow(dist,13))/dist;
        res[1] = -4*ewca*ndistv[1]*(6*pow((sigwca/pow(2,1.0/6.0)),6)/pow(dist,7) - 12*pow((sigwca/pow(2,1.0/6.0)),12)/pow(dist,13))/dist;
        res[2] = -4*ewca*ndistv[2]*(6*pow((sigwca/pow(2,1.0/6.0)),6)/pow(dist,7) - 12*pow((sigwca/pow(2,1.0/6.0)),12)/pow(dist,13))/dist;
    } else{
        res[0] = 0;
        res[1] = 0;
        res[2] = 0;
    }
}

void flj2p(double *x, double *y, double *z, double *rx, double *ry, double *rz){
    double comdif[3] = {-((x[0]+x[1])-(x[2]+x[3]))/2.0, -((y[0]+y[1])-(y[2]+y[3]))/2.0, -((z[0]+z[1])-(z[2]+z[3]))/2.0};
    double ncomdif[3] = {0, 0, 0};
    mic(comdif, ncomdif);
    double r = vec_size(ncomdif);
    if (r <= coffljc){
        double dudr = 4*eljc*(6*pow((sigljc/pow(2,1.0/6.0)),6)/pow(r,7) - 12*pow((sigljc/pow(2,1.0/6.0)),12)/pow(r,13))/r;
        rx[0] = -dudr/r/2*ncomdif[0]*(-1);
        ry[0] = -dudr/r/2*ncomdif[1]*(-1);
        rz[0] = -dudr/r/2*ncomdif[2]*(-1);
        rx[1] = -dudr/r/2*ncomdif[0]*(-1);
        ry[1] = -dudr/r/2*ncomdif[1]*(-1);
        rz[1] = -dudr/r/2*ncomdif[2]*(-1);
        rx[2] = -dudr/r/2*ncomdif[0];
        ry[2] = -dudr/r/2*ncomdif[1];
        rz[2] = -dudr/r/2*ncomdif[2];
        rx[3] = -dudr/r/2*ncomdif[0];
        ry[3] = -dudr/r/2*ncomdif[1];
        rz[3] = -dudr/r/2*ncomdif[2];
    } else{
        rx[0] = 0;
        ry[0] = 0;
        rz[0] = 0;
        rx[1] = 0;
        ry[1] = 0;
        rz[1] = 0;
        rx[2] = 0;
        ry[2] = 0;
        rz[2] = 0;
        rx[3] = 0;
        ry[3] = 0;
        rz[3] = 0;
    }
}

void fbonded(double *v1, double *v2, double *res){
    double distv[3] = {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
    double ndistv[3] = {0, 0, 0};
    mic(distv, ndistv);
    double dist = vec_size(ndistv);
    res[0] = -1*eb*(dist-r0b)/dist*distv[0];
    res[1] = -1*eb*(dist-r0b)/dist*distv[1];
    res[2] = -1*eb*(dist-r0b)/dist*distv[2];
}

void f_tot(double *x, double *y, double *z, double *Fx, double *Fy, double *Fz){
    if (US_on){
        // Umbrella sampling (2 MOLECULES ONLY)
        double f_U_dist_x[4];
        double f_U_dist_y[4];
        double f_U_dist_z[4];
        for (int ii=0; ii<4; ii++){
            f_U_dist_x[ii] = 0;
            f_U_dist_y[ii] = 0;
            f_U_dist_z[ii] = 0;
        }
        double fudx[4] = {x[0], x[1], x[2], x[3]};
        double fudy[4] = {y[0], y[1], y[2], y[3]};
        double fudz[4] = {z[0], z[1], z[2], z[3]};
        f_U_distance(fudx, fudy, fudz, f_U_dist_x, f_U_dist_y, f_U_dist_z, krU, r0U);


        // Umbrella sampling (2 MOLECULES ONLY)
        double f_U_q_x[4];
        double f_U_q_y[4];
        double f_U_q_z[4];
        for (int ii=0; ii<4; ii++){
            f_U_q_x[ii] = 0;
            f_U_q_y[ii] = 0;
            f_U_q_z[ii] = 0;
        }
        double fuqx[4] = {x[0], x[1], x[2], x[3]};
        double fuqy[4] = {y[0], y[1], y[2], y[3]};
        double fuqz[4] = {z[0], z[1], z[2], z[3]};
        f_U_q(fuqx, fuqy, fuqz, f_U_q_x, f_U_q_y, f_U_q_z, kqU, q0U);

        for (int ii=0; ii<4; ii++){
            Fx[ii] += f_U_dist_x[ii] + f_U_q_x[ii];
            Fy[ii] += f_U_dist_y[ii] + f_U_q_y[ii];
            Fz[ii] += f_U_dist_z[ii] + f_U_q_z[ii];
        }
    }
    // LJ, bonded, wca
    for (int i=0; i<2*n_mols; i++){
        if (i%4 == 0){
            double vb1[3] = {x[i], y[i], z[i]};
            double vb2[3] = {x[i+1], y[i+1], z[i+1]};
            double fbv[3] = {0, 0, 0};
            fbonded(vb1, vb2, fbv);
            //printf("bonded: %d\t%d\n", i, i+1);
            Fx[i] = Fx[i]+fbv[0];
            Fy[i] = Fy[i]+fbv[1];
            Fz[i] = Fz[i]+fbv[2];
            for (int j=0; j<2*n_mols; j++){
                if (j%4 ==2 || j%4 == 3){
                    double vlj1[3] = {x[i], y[i], z[i]};
                    double vlj2[3] = {x[j], y[j], z[j]};
                    double fljv[3] = {0, 0, 0};
                    flj(vlj1, vlj2, fljv);
                    //printf("lj: %d\t%d\n", i, j);
                    Fx[i] = Fx[i]+fljv[0];
                    Fy[i] = Fy[i]+fljv[1];
                    Fz[i] = Fz[i]+fljv[2];
                } else if(j != i && j!= i+1){
                    double vwca1[3] = {x[i], y[i], z[i]};
                    double vwca2[3] = {x[j], y[j], z[j]};
                    double fwcav[3] = {0, 0, 0};
                    fwca(vwca1, vwca2, fwcav);
                    //printf("wca: %d\t%d\n", i, j);
                    Fx[i] = Fx[i]+fwcav[0];
                    Fy[i] = Fy[i]+fwcav[1];
                    Fz[i] = Fz[i]+fwcav[2];
                }
            }
        }else if (i%4 == 1){
            double vb1[3] = {x[i], y[i], z[i]};
            double vb2[3] = {x[i-1], y[i-1], z[i-1]};
            double fbv[3]= {0 ,0, 0};
            fbonded(vb1, vb2, fbv);
            //printf("bonded: %d\t%d\n", i, i-1);
            Fx[i] = Fx[i]+fbv[0];
            Fy[i] = Fy[i]+fbv[1];
            Fz[i] = Fz[i]+fbv[2];
            for (int j=0; j<2*n_mols; j++){
                if (j%4 ==2 || j%4 == 3){
                    double vlj1[3] = {x[i], y[i], z[i]};
                    double vlj2[3] = {x[j], y[j], z[j]};
                    double fljv[3] = {0, 0, 0};
                    flj(vlj1, vlj2, fljv);
                    //printf("lj: %d\t%d\n", i, j);
                    Fx[i] = Fx[i]+fljv[0];
                    Fy[i] = Fy[i]+fljv[1];
                    Fz[i] = Fz[i]+fljv[2];
                } else if(j != i && j!= i-1){
                    double vwca1[3] = {x[i], y[i], z[i]};
                    double vwca2[3] = {x[j], y[j], z[j]};
                    double fwcav[3] = {0, 0, 0};
                    fwca(vwca1, vwca2, fwcav);
                    //printf("wca: %d\t%d\n", i, j);
                    Fx[i] = Fx[i]+fwcav[0];
                    Fy[i] = Fy[i]+fwcav[1];
                    Fz[i] = Fz[i]+fwcav[2];
                }
            }
        }else if(i%4 == 2){
            double vb1[3] = {x[i], y[i], z[i]};
            double vb2[3] = {x[i+1], y[i+1], z[i+1]};
            double fbv[3] = {0, 0, 0};
            fbonded(vb1, vb2, fbv);
            //printf("bonded: %d\t%d\n", i, i+1);
            Fx[i] = Fx[i]+fbv[0];
            Fy[i] = Fy[i]+fbv[1];
            Fz[i] = Fz[i]+fbv[2];
            for (int j=0; j<2*n_mols; j++){
                if (j%4 == 0 || j%4 == 1){
                    double vlj1[3] = {x[i], y[i], z[i]};
                    double vlj2[3] = {x[j], y[j], z[j]};
                    double fljv[3] = {0, 0, 0};
                    flj(vlj1, vlj2, fljv);
                    //printf("lj: %d\t%d\n", i, j);
                    Fx[i] = Fx[i]+fljv[0];
                    Fy[i] = Fy[i]+fljv[1];
                    Fz[i] = Fz[i]+fljv[2];
                } else if(j != i && j!= i+1){
                    double vwca1[3] = {x[i], y[i], z[i]};
                    double vwca2[3] = {x[j], y[j], z[j]};
                    double fwcav[3]= {0, 0, 0};
                    fwca(vwca1, vwca2, fwcav);
                    //printf("wca: %d\t%d\n", i, j);
                    Fx[i] = Fx[i]+fwcav[0];
                    Fy[i] = Fy[i]+fwcav[1];
                    Fz[i] = Fz[i]+fwcav[2];
                }
            }
        }else if (i%4 == 3){
            double vb1[3] = {x[i], y[i], z[i]};
            double vb2[3] = {x[i-1], y[i-1], z[i-1]};
            double fbv[3] = {0, 0, 0};
            fbonded(vb1, vb2, fbv);
            //printf("bonded: %d\t%d\n", i, i-1);
            Fx[i] = Fx[i]+fbv[0];
            Fy[i] = Fy[i]+fbv[1];
            Fz[i] = Fz[i]+fbv[2];
            for (int j=0; j<2*n_mols; j++){
                if (j%4 == 0 || j%4 == 1){
                    double vlj1[3] = {x[i], y[i], z[i]};
                    double vlj2[3] = {x[j], y[j], z[j]};
                    double fljv[3] = {0, 0, 0};
                    flj(vlj1, vlj2, fljv);
                    //printf("lj: %d\t%d\n", i, j);
                    Fx[i] = Fx[i]+fljv[0];
                    Fy[i] = Fy[i]+fljv[1];
                    Fz[i] = Fz[i]+fljv[2];
                } else if(j != i && j!= i-1){
                    double vwca1[3] = {x[i], y[i], z[i]};
                    double vwca2[3] = {x[j], y[j], z[j]};
                    double fwcav[3] = {0, 0, 0};
                    fwca(vwca1, vwca2, fwcav);
                    //printf("wca: %d\t%d\n", i, j);
                    Fx[i] = Fx[i]+fwcav[0];
                    Fy[i] = Fy[i]+fwcav[1];
                    Fz[i] = Fz[i]+fwcav[2];
                }
            }
        } else {
            printf("This should not have happened...");
        }
    }
    // Angle + LJc
    for (int i=0; i < n_mols; i++){
        for (int j = 0; j<n_mols; j++){
            if (i != j && (2*i)%4 != (2*j)%4){
                //printf("%d\t%d\n", i, j);
                double tempx[4] = {x[2*i], x[2*i+1], x[2*j], x[2*j+1]};
                double tempy[4] = {y[2*i], y[2*i+1], y[2*j], y[2*j+1]};
                double tempz[4] = {z[2*i], z[2*i+1], z[2*j], z[2*j+1]};
                double fanglevx[4] = {0, 0, 0, 0};
                double fanglevy[4] = {0, 0, 0, 0};
                double fanglevz[4] = {0, 0, 0, 0};
                double fljc2pvx[4] = {0, 0, 0, 0};
                double fljc2pvy[4] = {0, 0, 0, 0};
                double fljc2pvz[4] = {0, 0, 0, 0};
		// The Wall
                fangle(tempx, tempy, tempz, fanglevx, fanglevy, fanglevz, ea, r0a, sigra, q0a, sigqa);
                flj2p(tempx, tempy, tempz, fljc2pvx, fljc2pvy, fljc2pvz);
                //printf("angle: %d\t%d\n", 2*i, 2*j);
                // The Channel
                double fcvx[4] = {0, 0, 0, 0};
                double fcvy[4] = {0, 0, 0, 0};
                double fcvz[4] = {0, 0, 0, 0};
                fangle(tempx, tempy, tempz, fcvx, fcvy, fcvz, ec, r0c, sigrc, q0c, sigqc);
                // Tilting locked state
                double ftiltvxl[4] = {0, 0, 0, 0};
                double ftiltvyl[4] = {0, 0, 0, 0};
                double ftiltvzl[4] = {0, 0, 0, 0};
                fangle(tempx, tempy, tempz, ftiltvxl, ftiltvyl, ftiltvzl, etl, r0t, sigrt, q0tl, sigqt);
                Fx[2*i] = Fx[2*i] + fanglevx[0] + fljc2pvx[0] + fcvx[0]+ ftiltvxl[0];
                Fy[2*i] = Fy[2*i] + fanglevy[0] + fljc2pvy[0] + fcvy[0]+ ftiltvyl[0];
                Fz[2*i] = Fz[2*i] + fanglevz[0] + fljc2pvz[0] + fcvz[0]+ ftiltvzl[0];
                Fx[2*i+1] = Fx[2*i+1] + fanglevx[1] + fljc2pvx[1] + fcvx[1] + ftiltvxl[1];
                Fy[2*i+1] = Fy[2*i+1] + fanglevy[1] + fljc2pvy[1] + fcvy[1] + ftiltvyl[1];
                Fz[2*i+1] = Fz[2*i+1] + fanglevz[1] + fljc2pvz[1] + fcvz[1] + ftiltvzl[1];
            }
        }
    }
}

// ----------------------------------------------------------------------------------------------------------------------------------------------
// Integrator

void langevin_integrate(double *x, double *y, double *z, double tauq, gsl_rng *r){
    double forcex[2*n_mols];
    double forcey[2*n_mols];
    double forcez[2*n_mols];
    //forcex = (double*)malloc(2*n_mols*sizeof(double));
    //forcey = (double*)malloc(2*n_mols*sizeof(double));
    //forcez = (double*)malloc(2*n_mols*sizeof(double));
    for (int i=0; i<2*n_mols; i++){
        forcex[i] = 0;
        forcey[i] = 0;
        forcez[i] = 0;
    }
    f_tot(x, y, z, forcex, forcey, forcez);
    // for (int i=0; i<2*n_mols; i++){
    //     printf("%f\t", forcex[i]);
    //     printf("%f\t", forcey[i]);
    //     printf("%f\n", forcez[i]);
    // }
    for(int i = 0; i<2*n_mols; i++){
        double ksix = gsl_ran_gaussian(r,1);
        double ksiy = gsl_ran_gaussian(r,1);
        double ksiz = gsl_ran_gaussian(r,1);
        x[i] = x[i] + tauq*forcex[i] + sqrt(2*tauq)*ksix;
        y[i] = y[i] + tauq*forcey[i] + sqrt(2*tauq)*ksiy;
        z[i] = z[i] + tauq*forcez[i] + sqrt(2*tauq)*ksiz;
    }
    //free(forcex);
    //free(forcey);
    //free(forcez);
}
// ----------------------------------------------------------------------------------------------------------------------------------------------
// Input/Output

void write_output(char *filename, unsigned long long step, double *x, double *y, double *z, double *mass){
    double xn[2*n_mols];
    double yn[2*n_mols];
    double zn[2*n_mols];
    for (int i=0; i<2*n_mols; i++){
        xn[i] = pbc_cord(x[i]);
        yn[i] = pbc_cord(y[i]);
        zn[i] = pbc_cord(z[i]);
    }
    FILE *handle = fopen(filename, "a");
    fprintf(handle, "ITEM: TIMESTEP\n");
    fprintf(handle, "%lld\n", step);
    fprintf(handle, "ITEM: NUMBER OF ATOMS\n");
    fprintf(handle, "%d\n", 2*n_mols);
    fprintf(handle, "ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(handle, "%f\t%f\n", -box_size*rd/2,box_size*rd/2);
    fprintf(handle, "%f\t%f\n", -box_size*rd/2,box_size*rd/2);
    fprintf(handle, "%f\t%f\n", -box_size*rd/2,box_size*rd/2);
    fprintf(handle, "ITEM: ATOMS id mol type mass x y z\n");
    for(int i=0; i<2*n_mols; i++){
        fprintf(handle, "%d\t%d\t%d\t%f\t%f\t%f\t%f\n", i+1, (int) ceil((i+1)/2.0), (i%4)+1, mass[i], xn[i]*rd, yn[i]*rd, zn[i]*rd);
    }
    fclose(handle);
}

void write_colvars(char *filename, unsigned long long step, double rc, double q) {
    FILE *handle = fopen(filename, "a");
    fprintf(handle, "%lld\t%f\t%f\n", step, rc*rd, q);
    fclose(handle);
}

void write_colvars_L(char *filename, unsigned long long step, int L) {
    FILE *handle = fopen(filename, "a");
    fprintf(handle, "%lld\t%d\n", step, L);
    fclose(handle);
}

int read_restart(char *fname, double *x, double *y, double *z){
    FILE *handle = fopen(fname, "r");
    if(handle == NULL){
        printf("Error reading file!\n");
        return 1;
    }
    int inp_n_mol;
    fscanf(handle, "%d", &inp_n_mol);
    if(inp_n_mol != n_mols){
        printf("Wrong input file!\n");
        return 1;
    }
    for(int i=0; i<2*n_mols; i++){
        fscanf(handle, "%lf", &x[i]);
        }
    for(int i=0; i<2*n_mols; i++){
        fscanf(handle, "%lf", &y[i]);
        }
    for(int i=0; i<2*n_mols; i++){
        fscanf(handle, "%lf", &z[i]);
    }
    fclose(handle);
    return 0;
}
