#include "ns2d.hh"
#include <mpi.h>


int main(int argc, char *argv[]) {

    // use this flag to enable/disable density validation checks
    bool check_density = false;

   // Initialize MPI
    MPI_Init(&argc, &argv);

    Params params;
    parse_args(argc, argv, params);

    const int N = params.N;
    const double h = (2.0*M_PI) / N;
    const double nu = 1.0 / params.Re;
    const int max_size = (N+2)*(N+2);
    // Use this to control the number of warm-up steps to ignore in timing until the node is sufficiently warmed up and the timings are stable
    const int warmUp = 2; 

    // measure the wall clock time for the entire simulation
    // and the time for each step
    double t = 0.0;
    int step = 0;

    // accumulated phase timers (steps >= warmUp only)
    double t_halo = 0.0, t_rhs = 0.0, t_rk = 0.0, t_step = 0.0;
    int timed_steps = 0;
    double t0, t_step_start;

    
    printf("========================================\n");
    printf("  2D Compressible Isothermal NS\n");
    printf("========================================\n");
    printf("  N         = %d\n",   N);
    printf("  Re        = %.4f\n", params.Re);
    printf("  T         = %.4f\n", params.T);
    printf("  eps       = %.4f\n", params.eps);
    printf("  nu        = %.4f\n", nu);
    printf("  dt_factor = %.4f\n", params.dt_factor);
    printf("========================================\n");
    
    // Allocate memory for the fields and their derivatives
    double *rho = (double *)calloc(max_size, sizeof(double));
    double *u = (double *)calloc(max_size, sizeof(double));
    double *v = (double *)calloc(max_size, sizeof(double));

    double *rho1 = (double *)calloc(max_size, sizeof(double));
    double *u1 = (double *)calloc(max_size, sizeof(double));
    double *v1 = (double *)calloc(max_size, sizeof(double));

    double *rho2 = (double *)calloc(max_size, sizeof(double));
    double *u2 = (double *)calloc(max_size, sizeof(double));
    double *v2 = (double *)calloc(max_size, sizeof(double));

    double *rho_new = (double *)calloc(max_size, sizeof(double));
    double *u_new = (double *)calloc(max_size, sizeof(double));
    double *v_new = (double *)calloc(max_size, sizeof(double));

    double *rhs_rho = (double *)calloc(max_size, sizeof(double));
    double *rhs_u = (double *)calloc(max_size, sizeof(double));
    double *rhs_v = (double *)calloc(max_size, sizeof(double));

    // Initialize the fields with the Taylor-Green vortex
    initialize_TG(rho, u, v, h, N);
    double dt = 0.0;

    // Compute the initial time step based on CFL condition
    dt = compute_dt(u, v, N, h, nu, params.dt_factor);

    if (check_density) {
        if (!validate_density(rho, N, "initialization", step, 0)) {
            return EXIT_FAILURE;
        }
    }

    // Main time-stepping loop
    double t_total_start = MPI_Wtime();
    
    while (t < params.T) {

        if(step % 100 == 0) {
            double mesh_Re = get_mesh_Re(u, v, N, h, nu);
            double velocity_rms = get_velocity_rms(u, v, N);
            double density_rms = get_density_rms(rho, N);

            if(step == 0){
                printf("========================================\n");
                printf("  Step, Time, dt, Mesh Re, Velocity RMS, Density RMS\n");
                printf("========================================\n");
            }
            printf("%d, %.8f, %.8f, %.8f, %.8f, %.8f\n", step, t, dt, mesh_Re, velocity_rms, density_rms);
            fflush(stdout);
        }

        // do not step past the final time
        if (t + dt > params.T) dt = params.T - t;

        t_step_start = MPI_Wtime();

        // sub-step 1; fill halos, compute RHS, and update fields
        t0 = MPI_Wtime();
        fill_ghost_cells(rho, N);
        fill_ghost_cells(u, N);
        fill_ghost_cells(v, N);
        if (step >= warmUp) t_halo += MPI_Wtime() - t0;

        t0 = MPI_Wtime();
        compute_rhs(rho, u, v, rhs_rho, rhs_u, rhs_v, h, N, nu);
        if (step >= warmUp) t_rhs += MPI_Wtime() - t0;

        t0 = MPI_Wtime();
        rk3_step_1(rho, rho1, rhs_rho, h, N, dt);
        rk3_step_1(u, u1, rhs_u, h, N, dt);
        rk3_step_1(v, v1, rhs_v, h, N, dt);
        if (step >= warmUp) t_rk += MPI_Wtime() - t0;

        if (check_density) {
            if (!validate_density(rho1, N, "RK3 stage 1", step, 0)) {
                return EXIT_FAILURE;
            }
        }

        // sub-step 2; fill halos, compute RHS, and update fields
        t0 = MPI_Wtime();
        fill_ghost_cells(rho1, N);
        fill_ghost_cells(u1, N);
        fill_ghost_cells(v1, N);
        if (step >= warmUp) t_halo += MPI_Wtime() - t0;

        t0 = MPI_Wtime();
        compute_rhs(rho1, u1, v1, rhs_rho, rhs_u, rhs_v, h, N, nu);
        if (step >= warmUp) t_rhs += MPI_Wtime() - t0;

        t0 = MPI_Wtime();
        rk3_step_2(rho, rho1, rho2, rhs_rho, h, N, dt);
        rk3_step_2(u, u1, u2, rhs_u, h, N, dt);
        rk3_step_2(v, v1, v2, rhs_v, h, N, dt);
        if (step >= warmUp) t_rk += MPI_Wtime() - t0;

        if (check_density) {
            if (!validate_density(rho2, N, "RK3 stage 2", step, 0)) {
                return EXIT_FAILURE;
            }
        }

        // sub-step 3; fill halos, compute RHS, and update fields
        t0 = MPI_Wtime();
        fill_ghost_cells(rho2, N);
        fill_ghost_cells(u2, N);
        fill_ghost_cells(v2, N);
        if (step >= warmUp) t_halo += MPI_Wtime() - t0;

        t0 = MPI_Wtime();
        compute_rhs(rho2, u2, v2, rhs_rho, rhs_u, rhs_v, h, N, nu);
        if (step >= warmUp) t_rhs += MPI_Wtime() - t0;

        t0 = MPI_Wtime();
        rk3_step_3(rho, rho2, rho_new, rhs_rho, h, N, dt);
        rk3_step_3(u, u2, u_new, rhs_u, h, N, dt);
        rk3_step_3(v, v2, v_new, rhs_v, h, N, dt);
        if (step >= warmUp) t_rk += MPI_Wtime() - t0;

        if (check_density) {
            if (!validate_density(rho_new, N, "RK3 stage 3", step, 0)) {
                return EXIT_FAILURE;
            } 
        }

        // Update the fields for the next time step
        double *tmp;
        tmp = rho; rho = rho_new; rho_new = tmp;
        tmp = u;   u = u_new;     u_new = tmp;
        tmp = v;   v = v_new;     v_new = tmp;

        // Update the time and step
        t += dt;
        step++;

        // Compute the next time step based on CFL condition
        dt = compute_dt(u, v, N, h, nu, params.dt_factor);

        if (step - 1 >= warmUp) {
            t_step += MPI_Wtime() - t_step_start;
            timed_steps++;
        }
    }

    double t_total = MPI_Wtime() - t_total_start;
    double l2 = compute_l2_error(u, v, N, h, t, nu);

    printf("========================================\n");
    printf("  Simulation completed\n");
    printf("  Total steps: %d\n", step);
    printf("  Total simulation time: %.6f seconds\n", t);
    printf("  Wall clock time: %.6f seconds\n", t_total);
    printf("  L2 error: %.6e\n", l2);
    printf("========================================\n");

    if (timed_steps > 0) {
        double n = (double)timed_steps;
        printf("  Average time per step over %d steps (first %d skipped)\n", timed_steps, warmUp);
        printf("  Halo exchange : %12.3f us  (%5.1f %%)\n", 1e6 * t_halo / n, 100.0 * t_halo / t_step);
        printf("  RHS           : %12.3f us  (%5.1f %%)\n", 1e6 * t_rhs / n, 100.0 * t_rhs / t_step);
        printf("  RK3 update    : %12.3f us  (%5.1f %%)\n", 1e6 * t_rk / n, 100.0 * t_rk / t_step);
        printf("  Full step     : %12.3f us\n", 1e6 * t_step / n);
        printf("========================================\n");
        printf("TIMING,N=%d,halo_us=%.4f,rhs_us=%.4f,rk_us=%.4f,step_us=%.4f,wall_s=%.6f,l2=%.6e\n",
               N, 1e6 * t_halo / n, 1e6 * t_rhs / n, 1e6 * t_rk / n, 1e6 * t_step / n, t_total, l2);
    }

    free(rho);
    free(u);
    free(v);
    free(rho1);
    free(u1);
    free(v1);
    free(rho2);
    free(u2);
    free(v2);
    free(rho_new);
    free(u_new);
    free(v_new);  
    free(rhs_rho);
    free(rhs_u);
    free(rhs_v);

    MPI_Finalize();
    return 0;

}