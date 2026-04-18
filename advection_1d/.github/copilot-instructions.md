# Copilot instructions for `advection_1d`

## Build, test, and lint commands

Run commands from the `advection_1d/` directory.

```bash
# Configure (Debug; tests enabled by default)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug

# Build solver and tests
cmake --build build --parallel

# Run full test suite
ctest --test-dir build --output-on-failure

# Run one CTest target
ctest --test-dir build -R test_operator --output-on-failure

# Run one GoogleTest case directly
./build/tests/test_reference_element --gtest_filter="*DMatrix*"
```

No dedicated lint target is defined in CMake or project docs.

## High-level architecture

- The solver is organized as a templated DG pipeline in `headers/`, wired in `src/main.cpp`: **Mesh** → **ReferenceElement** → **NumericalFlux** → **Operator (RHS)** → **Integrator** → **VTK output**.
- `Mesh.hpp` wraps a 1D TNL mesh and owns DG-specific geometry/adjacency helpers (element sizes, face neighbors, normals, boundary detection), with factories `uniform(...)` and `readVTK(...)`.
- `ReferenceElement.hpp` builds nodal DG operators (GLL nodes/weights, Vandermonde, `Dr`, `LIFT`) using TNL matrices plus Eigen-based inversion helpers.
- `Operator.hpp` computes the semi-discrete RHS per element from volume + interface terms and delegates interface treatment to `NumericalFlux` implementations.
- `Integrator.hpp` provides multiple explicit RK variants (`ERK`, `LSERK`, `SSPRK`) that advance `FieldVector` using a shared RHS callback signature.
- `IO.hpp` writes element-local DG states to legacy ASCII VTK (`*_000000.vtk` style time series).

## Key codebase conventions

- Keep new numerical code in the `DG` namespace and follow the existing template defaults (`Real=double`, `Device=TNL::Devices::Host`, `Index=int`).
- `FieldVector` storage is contiguous and element-major: index as `i + k * Np` (or use `elementPtr(k)` for element-local access).
- Face indexing is vertex-based in 1D (`0..K`), with `Mesh::BOUNDARY_FACE == -1` and neighbor queries done through `leftCellOfFace/rightCellOfFace`.
- Boundary behavior in `Operator::computeRHS` is currently periodic (left boundary reads from last element, right boundary reads from first). Preserve this unless explicitly changing boundary conditions.
- Time integrators expect RHS callbacks with signature `void(const Field&, Field&, const Real&)`; reuse `Operator::rhsFunction()` instead of creating ad-hoc adapters.
- When adding tests, register each new test executable in `tests/CMakeLists.txt` through `dg_add_test(...)` so `gtest_discover_tests` exposes it to `ctest`.
- `README.md` is the canonical project usage doc; `.github/copilot-code-review-instructions.md` exists separately for review-focused guidance on `*.cpp`/`*.hpp` files.
