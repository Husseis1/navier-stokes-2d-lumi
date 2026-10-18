#ifndef NS2D_HH
#define NS2D_HH

#include <math.h>
#include <cmath>
#include <vector>
#include <cfloat>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>

#define IDX(i, j, N) ((j+1)*(N+2) + (i+1))


typedef struct {
    int N;
    double Re;
    double T;
    double eps;
    double dt_factor;
} Params;

/*Utility functions*/
void parse_args(int argc, char *argv[], Params &params);
bool validate_density(const double *rho, int N, const char *stage, int step, int rank);
double get_mesh_Re(double *u, double *v, int N, double h, double nu);
double get_velocity_rms(double *u, double *v, int N);
double get_density_rms(double *rho, int N);
double compute_l2_error(double *u, double *v, int N, double h, double t, double nu);


/* TODO: Implement the following functions */

void rk3_step_1(double *f, double *f1, double *rhs, double h, int N, double dt);
void rk3_step_2(double *f, double *f1, double *f2, double *rhs, double h, int N, double dt);
void rk3_step_3(double *f, double *f2, double *f_new, double *rhs, double h, int N, double dt);

double compute_dt(const double *u, const double *v, int N, double h, double nu, double dt_factor);

void compute_rhs(double *rho, double *u, double *v, double *rhs_rho, double *rhs_u, double *rhs_v, double h, int N, double nu);

void initialize_TG(double *rho, double *u, double *v, double h, int N);

void fill_ghost_cells(double *f, int N);




#endif