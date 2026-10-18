#include "ns2d.hh"

/*Utility functions*/
void parse_args(int argc, char *argv[], Params &params) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--N" && i + 1 < argc) {
            params.N = std::stoi(argv[++i]);
        } else if (arg == "--Re" && i + 1 < argc) {
            params.Re = std::stod(argv[++i]);
        } else if (arg == "--T" && i + 1 < argc) {
            params.T = std::stod(argv[++i]);
        } else if (arg == "--eps" && i + 1 < argc) {
            params.eps = std::stod(argv[++i]);
        } else if (arg == "--dt_factor" && i + 1 < argc) {
            params.dt_factor = std::stod(argv[++i]);
        } else {
            fprintf(stderr, "Unknown or incomplete argument: %s\n", arg.c_str());
            fprintf(stderr, "Usage: %s --N <grid_size> --Re <Reynolds_number> --T <final_time> --eps <tolerance> --dt_factor <time_step_factor>\n", argv[0]);
            fflush(stderr);
            exit(EXIT_FAILURE);
        }
    }
}

bool validate_density(const double *rho, int N, const char *stage, int step, int rank) {
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            const double value = rho[IDX(i, j, N)];
            if (!std::isfinite(value) || value <= 0.0) {
                if (rank == 0) {
                    fprintf(stderr,
                            "Invalid density at step %d during %s: rho[%d,%d]=%.17e\n",
                            step, stage, i, j, value);
                    fflush(stderr);
                }
                return false;
            }
        }
    }

    return true;
}

double get_mesh_Re(double *u, double *v, int N, double h, double nu) {
    double max_magnitude = 0.0;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            double magnitude = u[IDX(i, j, N)] * u[IDX(i, j, N)] + v[IDX(i, j, N)] * v[IDX(i, j, N)];
            if (magnitude > max_magnitude) {
                max_magnitude = magnitude;
            }
        }
    }

    // Compute the maximum velocity magnitude
    double max_velocity = std::sqrt(max_magnitude);

    // Compute the mesh Reynolds number
    return (max_velocity * h) / nu;
}

double get_velocity_rms(double *u, double *v, int N) {
    double sum_sq = 0.0;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            sum_sq += u[IDX(i, j, N)] * u[IDX(i, j, N)] + v[IDX(i, j, N)] * v[IDX(i, j, N)];
        }
    }
    return std::sqrt(sum_sq / (double)(N * N));
}

double get_density_rms(double *rho, int N) {
    double sum_sq = 0.0;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            sum_sq += rho[IDX(i, j, N)] * rho[IDX(i, j, N)];
        }
    }
    return std::sqrt(sum_sq / (double)(N * N));
}

// Computes the L2 error between the numerical solution and the exact solution for the Taylor-Green vortex problem
double compute_l2_error(double *u, double *v, int N, double h, double t, double nu) {
    double l2_error = 0.0;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            double x = (i + 0.5) * h;
            double y = (j + 0.5) * h;
            double exact_u = sin(x) * cos(y) * exp(-2.0 * nu * t);
            double exact_v = -cos(x) * sin(y) * exp(-2.0 * nu * t);
            double error_u = u[IDX(i, j, N)] - exact_u;
            double error_v = v[IDX(i, j, N)] - exact_v;
            l2_error += ((error_u * error_u) + (error_v * error_v));
        }
    }
    l2_error = sqrt(l2_error / (N * N));
    return l2_error;
}



void rk3_step_1(double *f, double *f1, double *rhs, double h, int N, double dt){
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            const int k = IDX(i, j, N);
            f1[k] = f[k] + dt * rhs[k];
        }
    }
}
void rk3_step_2(double *f, double *f1, double *f2, double *rhs, double h, int N, double dt){
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            const int k = IDX(i, j, N);
            f2[k] = 0.75 * f[k] + 0.25 * f1[k] + 0.25 * dt * rhs[k];
        }
    }
}
void rk3_step_3(double *f, double *f2, double *f_new, double *rhs, double h, int N, double dt){
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            const int k = IDX(i, j, N);
            f_new[k] = (1.0 / 3.0) * f[k] + (2.0 / 3.0) * f2[k] + (2.0 / 3.0) * dt * rhs[k];
        }
    }
}

double compute_dt(const double *u, const double *v, int N, double h, double nu, double dt_factor){
    const double cs = 1.0;
    double umax = 0.0, vmax = 0.0;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            umax = fmax(umax, fabs(u[IDX(i, j, N)]));
            vmax = fmax(vmax, fabs(v[IDX(i, j, N)]));
        }
    }
    return dt_factor * h * h / (4.0 * nu + h * (umax + vmax + 2.0 * cs));
}


void initialize_TG(double *rho, double *u, double *v, double h, int N) {
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            double x = (i + 0.5) * h;
            double y = (j + 0.5) * h;
            rho[IDX(i, j, N)] = 1.0;
            u[IDX(i, j, N)] = sin(x) * cos(y);
            v[IDX(i, j, N)] = -cos(x) * sin(y);
        }
    }
}


// Computes the right-hand side of the vorticity equation using second order finite differences
void compute_rhs(double *rho, double *u, double *v, double *rhs_rho, double *rhs_u, double *rhs_v, double h, int N, double nu) {
    const double cs2 = 1.0;
    const double inv2h = 1.0 / (2.0 * h);
    const double invh2 = 1.0 / (h * h);
    const double inv4h2 = 1.0 / (4.0 * h * h);

    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            const int c  = IDX(i,     j,     N);
            const int e  = IDX(i + 1, j,     N);
            const int w  = IDX(i - 1, j,     N);
            const int n  = IDX(i,     j + 1, N);
            const int s  = IDX(i,     j - 1, N);
            const int ne = IDX(i + 1, j + 1, N);
            const int se = IDX(i + 1, j - 1, N);
            const int nw = IDX(i - 1, j + 1, N);
            const int sw = IDX(i - 1, j - 1, N);

            const double r = rho[c];
            const double uc = u[c];
            const double vc = v[c];

            const double drdx = (rho[e] - rho[w]) * inv2h;
            const double drdy = (rho[n] - rho[s]) * inv2h;
            const double dudx = (u[e] - u[w]) * inv2h;
            const double dudy = (u[n] - u[s]) * inv2h;
            const double dvdx = (v[e] - v[w]) * inv2h;
            const double dvdy = (v[n] - v[s]) * inv2h;

            const double d2udx2 = (u[e] - 2.0 * uc + u[w]) * invh2;
            const double d2udy2 = (u[n] - 2.0 * uc + u[s]) * invh2;
            const double d2vdx2 = (v[e] - 2.0 * vc + v[w]) * invh2;
            const double d2vdy2 = (v[n] - 2.0 * vc + v[s]) * invh2;
            const double d2udxdy = (u[ne] - u[se] - u[nw] + u[sw]) * inv4h2;
            const double d2vdxdy = (v[ne] - v[se] - v[nw] + v[sw]) * inv4h2;

            const double Sxx = (2.0 / 3.0) * dudx - (1.0 / 3.0) * dvdy;
            const double Syy = (2.0 / 3.0) * dvdy - (1.0 / 3.0) * dudx;
            const double Sxy = 0.5 * (dudy + dvdx);

            rhs_rho[c] = -uc * drdx - vc * drdy - r * (dudx + dvdy);

            rhs_u[c] = -uc * dudx - vc * dudy
                       - (cs2 / r) * drdx
                       + nu * (d2udx2 + d2udy2)
                       + (nu / 3.0) * (d2udx2 + d2vdxdy)
                       + (2.0 * nu / r) * (Sxx * drdx + Sxy * drdy);

            rhs_v[c] = -uc * dvdx - vc * dvdy
                       - (cs2 / r) * drdy
                       + nu * (d2vdx2 + d2vdy2)
                       + (nu / 3.0) * (d2vdy2 + d2udxdy)
                       + (2.0 * nu / r) * (Sxy * drdx + Syy * drdy);
        }
    }
}

// Fill ghost cells for periodic boundary conditions
void fill_ghost_cells(double *f, int N) {
    for (int j = 0; j < N; ++j) {
        f[IDX(-1, j, N)] = f[IDX(N - 1, j, N)];
        f[IDX(N, j, N)] = f[IDX(0, j, N)];
    }
    for (int i = -1; i <= N; ++i) {
        f[IDX(i, -1, N)] = f[IDX(i, N - 1, N)];
        f[IDX(i, N, N)] = f[IDX(i, 0, N)];
    }
}