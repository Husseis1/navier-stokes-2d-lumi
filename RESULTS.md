# Part 1 – Serial baseline: results and discussion

## Setup
- LUMI-C, 1 MPI rank (AMD EPYC 7763 "Milan"), Cray clang 19.0.0 (PrgEnv-cray 8.6.0), `-O3`
- Taylor–Green vortex, Re = 10 (ν = 0.1), T = 10, dt_factor = 0.4
- Timers: `MPI_Wtime` around each phase. The first 2 steps are excluded as warm-up.
  Every value is the median of 3 runs. Run-to-run variation was below ~5 %.

Run command:
```bash
srun --export=ALL build/part1_ns --N <N> --Re 10 --T 10.0 --eps 1e-6 --dt_factor 0.4
```

## Results

Average time per RK3 time step (µs), with each phase's share of the full step:

| N | steps | halo exchange | RHS | RK3 update | other (dt) | full step | wall time [s] | L2 error |
|---|---|---|---|---|---|---|---|---|
| 32  | 619  | 0.41 (3.3 %) | 7.96 (65.4 %)   | 2.90 (23.8 %)  | 0.91 (7.4 %)  | 12.17  | 0.0076 | 1.33e-2 |
| 64  | 1755 | 0.69 (1.6 %) | 29.19 (66.5 %)  | 11.02 (25.1 %) | 3.00 (6.8 %)  | 43.90  | 0.077  | 1.40e-2 |
| 128 | 5585 | 2.05 (1.0 %) | 128.63 (64.8 %) | 57.02 (28.7 %) | 10.75 (5.4 %) | 198.44 | 1.11   | 1.41e-2 |

"Other" is the full step minus the three timed phases. It is mostly the CFL `dt` computation, a max-reduction over the grid.

Growth factor when N doubles:

| | 32 → 64 | 64 → 128 | ideal |
|---|---|---|---|
| RHS             | 3.67× | 4.41× | 4× (O(N²)) |
| RK3 update      | 3.81× | 5.17× | 4× (O(N²)) |
| halo exchange   | 1.70× | 2.96× | 2× (O(N))  |
| full step       | 3.61× | 4.52× | 4×         |
| number of steps | 2.84× | 3.18× | 2–4×       |
| total wall time | 10.2× | 14.4× | 8–16×      |

## Discussion

**Which phase dominates?**
The RHS evaluation dominates at every N, taking about 65 % of each time step. It is the only phase with substantial arithmetic. For every grid point it evaluates a 9-point stencil for three fields: six first derivatives, four second derivatives, two mixed derivatives, the strain tensor, and the three right-hand sides. That is roughly 100 floating-point operations per point, including a division by ρ. The cost is about 7–8 ns per grid point per step, almost independent of N.

The RK3 update is second, at 24–29 %. It is a pure streaming kernel: each point reads 2–3 values, does about 3 flops, and writes one value. It is therefore limited by memory bandwidth rather than computation. The halo exchange costs 1–3 %, and the `dt` reduction 5–7 %.

**Scaling of each phase with N.**
- **RHS and RK3** loop over all N² interior points, so their cost should grow 4× per doubling of N. The RHS follows this closely (3.7×, 4.4×).
- **The RK3 update grows faster at large N** (5.2× from 64 to 128), and its share rises from 24 % to 29 %. The 15 field arrays need about 0.5 MB at N = 64, which fits in Milan's 512 KB L2 cache. At N = 128 they need about 2 MB, which spills into the L3 cache. Higher-latency, lower-bandwidth memory hurts the memory-bound RK3 update the most. The compute-heavy RHS is less sensitive, so its share drops slightly.
- **The halo exchange** touches only the 4(N+1) boundary cells per field, so it scales as O(N): 1.7× from 32 to 64. From 64 to 128 it grows 3×. The left and right ghost columns are strided accesses (stride N+2, so each copied value sits on a different cache line), and these become more expensive as the arrays leave the L2 cache. Because it is O(N) against O(N²) for everything else, its share keeps shrinking: 3.3 % → 1.6 % → 1.0 %.
- **The full step** therefore scales close to O(N²): 3.6×, then 4.5×.

**Cost of the full simulation.**
The time step is set by the CFL condition, `dt = C h² / (4ν + h(|u|max + |v|max + 2c_s))`. When the advective term dominates, dt ∝ h and the step count grows 2× per doubling of N. When the viscous term dominates, dt ∝ h² and the step count grows 4×. At Re = 10 and these resolutions the two terms are similar in size, so the step count grows 2.8× and then 3.2×, trending towards 4×. Total wall time therefore grows between O(N³) and O(N⁴): measured 10.2× and 14.4× per doubling. Going from N = 32 to N = 128 costs about 150× more wall time for a 16× larger grid. At higher resolution or lower Re the diffusive limit takes over and the cost approaches O(N⁴). This is the motivation for the parallelisation in the following parts.

**Outlook for Parts 2–4.**
In serial the halo exchange is a cheap periodic copy. With MPI domain decomposition it becomes real communication with latency and bandwidth costs, while each rank's share of the O(N²) work shrinks. The halo's share of the step will therefore grow with the number of ranks, which is the main factor limiting strong scaling.

The RHS and RK3 kernels are natural targets for optimisation. Fusing the RK3 update into the RHS loop would save one full read/write sweep per stage, which is valuable because the RK3 update is memory-bound. All phases map directly onto one-thread-per-point GPU kernels.

## Verification

**Full amplitude.** At T = 10 the L2 error against the exact solution is about 1.3–1.4·10⁻² and does not decrease with N. At T = 1 it is about 0.135 for every N. This is not a discretisation error:
- The initial condition has |u|max = 1 = c_s, so the flow is at Mach ~1.
- The uniform initial density does not supply the pressure field the incompressible vortex needs.
- The flow therefore launches acoustic waves, and the density RMS drifts away from 1 (1.027 at t = 1). The incompressible exact solution only holds approximately (ρ ≈ 1).

**Small amplitude.** To verify the discretisation itself, the same solver was run with the initial velocity scaled by A = 10⁻⁶, where compressible and nonlinear effects are negligible. The L2 error relative to the exact solution was:

| N | 16 | 32 | 64 | 128 |
|---|---|---|---|---|
| rel. L2 error | 7.56e-4 | 1.87e-4 | 4.65e-5 | 1.16e-5 |
| ratio | – | 4.05 | 4.01 | 4.00 |

The error falls by 4× per halving of h, which confirms second-order spatial accuracy.

## Notes on the handout
1. **CFL formula.** Simplifying the second line gives `dt = C h² / (4ν + h(|u|max + |v|max + 2c_s))`. The handout's last line leaves out the factor h in front of 2c_s, which makes it dimensionally inconsistent. The code uses the corrected form.
2. **y-momentum, ν/3 term.** The handout writes ν/3 (∂²v/∂x² + ∂²u/∂x∂y). The physically correct term is the y-component of ν/3 ∇(∇·u), which is ν/3 (∂²u/∂x∂y + ∂²v/∂y²). The code uses ∂²v/∂y². For the Taylor–Green field both versions give the same result.
3. **LUMI environment.** The executable needed `export LD_LIBRARY_PATH=/opt/cray/libfabric/1.22.0/lib64:$LD_LIBRARY_PATH` and `srun --export=ALL`, otherwise it fails at start-up with `libfabric.so.1: cannot open shared object file`.
