# 2D Compressible Navier–Stokes Solver (C++ / MPI, LUMI)

A finite-difference solver for the 2D compressible (weakly compressible, isothermal) Navier–Stokes equations on a periodic grid. Written in C++ with MPI. Benchmarked on the **LUMI supercomputer** (LUMI-C, AMD EPYC 7763).

This repository holds the **serial baseline** (Part 1) of a four-part project that takes the solver from a serial code to a distributed GPU implementation.

## What it does

- **Test case:** Taylor–Green vortex, with error measured against the analytical solution
- **Discretisation:** second-order central differences (9-point stencil), ghost-cell periodic boundaries
- **Time integration:** 3-stage Runge–Kutta (RK3) with an adaptive CFL time step
- **Instrumentation:** `MPI_Wtime` timers for each phase: halo exchange, RHS evaluation, RK3 update and the full time step

## Key results

Full write-up: [`RESULTS.md`](RESULTS.md)

- The **RHS evaluation dominates** at about 65 % of each step and costs about 7–8 ns per grid point. The **RK3 update** takes 24–29 %. It is memory-bound, and its share grows once the working set spills out of L2 cache (N = 128).
- Each step scales close to **O(N²)**. Total wall time scales between **O(N³) and O(N⁴)** because the CFL limit shifts from advective to diffusive.
- **Second-order accuracy verified.** In a small-amplitude run, the L2 error falls 4.0× each time h is halved (N = 16 → 128).
- The analysis also **found and corrected two errors in the provided equations**: a dimensionally inconsistent CFL formula and a wrong derivative in the y-momentum viscous term.

| N | full step [µs] | RHS share | wall time [s] |
|---|---|---|---|
| 32  | 12.2  | 65 % | 0.008 |
| 64  | 43.9  | 67 % | 0.077 |
| 128 | 198.4 | 65 % | 1.11  |

## Build and run

```bash
mkdir build && cd build
cmake .. && make -j
srun build/part1_ns --N 64 --Re 10 --T 1.0 --eps 1e-6 --dt_factor 0.4
```

| Argument | Meaning |
|---|---|
| `--N` | grid points per direction |
| `--Re` | Reynolds number |
| `--T` | final time |
| `--eps` | tolerance |
| `--dt_factor` | CFL safety factor |

`example-cpu-job.sh` is a SLURM job script for LUMI.

## Roadmap

1. Serial baseline and performance analysis *(this repo)*
2. Parts 2–4: MPI-parallel and GPU versions *(in progress)*
