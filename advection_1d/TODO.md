# Solver manager implementation steps

This file describes the execution order for introducing a single management `Solver` class that ties all existing solver parts together.

## Step-by-step plan

1. Define the new solver API in `include/Solver.hpp`.
   - Add `SolverConfig` with run parameters (`K`, `N`, domain, `Tf`, `CFL`, output names/cadence, Padé options).
   - Keep numerical-method choices compile-time (no runtime enum switching).
   - Declare `DG::Solver<Real, Device, Index, FluxPolicy, IntegratorPolicy, ICPolicy>` public methods (`initialize()`, `computeStableDt()`, `run()`).

2. Introduce compile-time policy selection in build configuration.
   - Add CMake options for method selection (e.g., flux type, integrator type, initial-condition preset).
   - Map each option to compile definitions or alias types used by `main.cpp` and/or `Solver.hpp`.
   - Fail configuration on invalid option values.

3. Implement internal assembly in `Solver` using compile-time policies.
   - Build mesh and reference element from config.
   - Build physical/advection lambdas.
   - Instantiate numerical flux from `FluxPolicy`.
   - Build `Operator` with selected flux and physical flux.
   - Instantiate integrator from `IntegratorPolicy` using `op.rhsFunction()`.

4. Move initialization setup into solver-owned routines.
   - Convert current `main.cpp` initial-condition kernels into reusable policy-driven setup methods.
   - Initialize `FieldVector` through selected `ICPolicy`.
   - Keep existing defaults/behavior consistent with the current executable.

5. Centralize time-step computation in solver.
   - Compute `x_min` from GLL spacing and mesh jacobians.
   - Compute max wave speed from initialized state.
   - Compute stable `dt` via selected integrator policy’s `computeDt(...)`.

6. Implement solver-managed run and output orchestration.
   - Write initial frame (`writeTimeSeriesVTK`).
   - Run timestep loop until `Tf`.
   - Write final frame and status output.
   - If enabled in config, run Padé reconstruction and write Padé output.

7. Refactor `src/main.cpp` into a thin entrypoint.
   - Replace manual pipeline wiring with config creation and selected compile-time aliases.
   - Construct `DG::Solver<...>` and call `run()`.
   - Preserve current defaults and output naming.

8. Add solver-level tests.
   - Add new tests for solver construction and short smoke run.
   - Add compile-time selection coverage for at least one non-default flux/integrator combination.
   - Register tests with `dg_add_test(...)` in `tests/CMakeLists.txt`.

9. Update usage documentation.
   - Update `README.md` with manager-based flow and compile-time selection instructions via CMake options.

10. Final validation pass.
    - Build project and run full tests.
    - Ensure behavior remains consistent with current baseline unless intentionally changed.

## Constraints to preserve
- Keep `Operator`, numerical flux formulas, and integrator algorithms unchanged.
- Keep periodic boundary behavior as currently implemented.
- Keep output compatibility (`output_XXXXXX.vtk` and optional Padé outputs).
