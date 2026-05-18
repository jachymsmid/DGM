# DGM solver for 1D advection equation

A templated Discontinuous Galerkin solver for the 1D advection equation, implemented following *Nodal Discontinuous Galerkin Methods* by Hesthaven and Warburton. The solver uses TNL (Template Numerical Library) for mesh and vector operations and supports multiple explicit time integrators (ERK, LSERK, SSPRK) with various numerical flux implementations.

## Architecture

The solver is organized as a templated DG pipeline:

**Mesh** → **ReferenceElement** → **NumericalFlux** → **Operator (RHS)** → **Integrator** → **VTK output**

All components are templated with:
- `Real`: floating-point type (default: `double`)
- `Device`: TNL device (default: `Host`)
- `Index`: integer type for indexing (default: `int`)

## MeshConfig.hpp
This is a config struct for the TNL::Meshes::Mesh constructor. It inherits from the TNL::Meshes::Mesh::DefaultConfig, that has these templates:
- Topology = TNL::Meshes::Topologies::Edge
- Dimension = 1
- GlobalIndex = int
- LocalIndex = short

We also set a alias for the mesh
```cpp
using TNLMesh1D = TNL::Meshes::Mesh<Mesh1DConfig<Real, Index, short>, Device>;
```
Where `Real`, `Index` and `Device` are templated.

## Mesh.hpp

The `Mesh` class wraps a TNL 1D mesh and provides DG-specific geometric utilities:

**Members:**
- `mesh_`: underlying TNL unstructured mesh
- `numElements_`: number of elements (K)
- Additional geometry and adjacency helpers (element sizes, face neighbors, normals)

**Constructors:**
- Construct from an existing TNL mesh
- Factory `uniform(a, b, K)`: creates a uniform mesh on interval [a,b] with K elements
- Factory `readVTK(filename)`: reads mesh from a VTK file

**Key methods:**
- `numElements()`: returns K (number of elements)
- `numFaces()`: returns K+1 (number of faces in 1D)
- `leftVertex(k)`, `rightVertex(k)`: boundary coordinates of element k
- `elementSize(k)`, `jacobian(k)`: element metrics
- `leftCellOfFace(f)`, `rightCellOfFace(f)`: neighbor queries (returns -1 at domain boundary)
- `isBoundaryFace(f)`: test if face is at domain boundary
- `faceCoord(f)`: physical coordinate of face f
- `leftNormal()`, `rightNormal()`: return -1 and +1 respectively (1D normals)
- `tnlMesh()`: access underlying TNL mesh

## ReferenceElement.hpp


Builds and manages nodal DG operators for the reference element [-1, 1].

**Members:**
- `N_`: polynomial order
- `Np_`: number of degrees of freedom (= N + 1)
- `r_`: Gauss-Legendre-Lobatto (GLL) node positions
- `w_`: GLL weights
- `V_`: Vandermonde matrix (evaluates Legendre basis at GLL nodes)
- `Dr_`: derivative matrix (computes derivatives of nodal expansion)
- `LIFT_`: lifting matrix (projects surface fluxes into volume RHS)

**Key methods:**
- `order()`: returns polynomial degree N
- `numDOF()`: returns Np (number of nodes per element)
- `nodes()`, `weights()`: GLL node/weight access
- `V()`, `Dr()`, `LIFT()`: operator matrices
- `legendreP(n, x)`: evaluates Legendre polynomial P_n(x)
- `legendrePderiv(n, x)`: evaluates dP_n/dx

**Implementation notes:** Uses Eigen-based matrix inversion and TNL matrices for storage. GLL nodes are computed via root-finding on derivatives of Legendre polynomials.


## FieldVector.hpp

Element-major contiguous storage for DG solution fields.

**Members:**
- `K_`: number of elements
- `Np_`: degrees of freedom per element
- `data_`: TNL vector holding K*Np values

**Storage layout:** Values are ordered element-major: `data[i + k*Np]` is the i-th DOF of element k.

**Constructors:**
- Default: empty (K=0, Np=0)
- `FieldVector(K, Np)`: allocate for K elements with Np DOFs each, initialize to zero

**Key methods:**
- `numElements()`, `numDOF()`, `totalSize()`: query dimensions
- `elementPtr(k)`: host pointer to start of element k's DOFs
- `data()`, `data() const`: access underlying storage
- `copyFrom(other)`: deep copy from another FieldVector

**Element-local access:** Use `elementPtr(k)` to obtain a local pointer for element-wise operations, or index the storage directly as `data()[i + k*Np]`.


## NumericalFlux.hpp

Abstract interface and implementations for computing numerical fluxes at element interfaces.

**Interface:** 
```cpp
template<class Real>
struct NumericalFlux {
  virtual Real compute(Real u_minus, Real u_plus, Real n_outward) const = 0;
};
```

Returns the scalar numerical flux f* to use in the DG surface term, given:
- `u_minus`: interior (left) state
- `u_plus`: exterior (right) state  
- `n_outward`: outward normal from interior element (-1 for left face, +1 for right face)

**Implementations:**

1. **UpwindFlux**: Selects the upwind state based on advection speed sign
   - If `a(u) * n_outward >= 0`: use interior flux `f(u_minus)`
   - Otherwise: use exterior flux `f(u_plus)`

2. **LaxFriedrichsFlux** (Rusanov): Arithmetic flux plus jump penalty
   - `f* = 0.5(f(u−) + f(u+)) + 0.5·C·n_outward·(u− − u+)`
   - where `C = max(|a(u−)|, |a(u+)|)`

3. **GodunovFlux**: Exact solver-based flux (for Burgers' equation)
   - Solves the local Riemann problem exactly

4. **RoeFlux**: Linearized Riemann solver using Roe average
   - Uses intermediate state weighted by characteristic speeds

All fluxes are constructed with function objects for advection speed and physical flux evaluation.

## Operator.hpp

Computes the semi-discrete RHS for the DG discretization via volume and interface terms.

**Members:**
- `mesh_`: DG mesh (geometry and adjacency)
- `ref_`: reference element (nodal operators)
- `flux_`: numerical flux implementation
- `physFlux_`: physical flux function callback

**Key method:**
- `computeRHS(const FieldVector& u, FieldVector& rhs)`: 
  - Computes volume term: `−Dr · u` on each element
  - Computes interface term: lift of numerical flux jumps
  - Stores semi-discrete RHS in `rhs`

**RHS callback:** Provides `rhsFunction()` signature for time integrators:
```cpp
void rhsFunction(const FieldVector& u, FieldVector& rhs, const Real& t)
```

**Boundary handling:** Currently uses periodic boundary conditions. For other conditions, modify the neighbor lookups in `computeRHS()`.


## Integrator.hpp

Multiple explicit Runge–Kutta time integration schemes with a shared callback signature.

**RHS callback signature:**
```cpp
void(const FieldVector& u, FieldVector& rhs, const Real& t)
```

**Integrator classes:**

1. **ERK** (Explicit Runge–Kutta)
   - Standard p-stage explicit RK of order p (p ∈ {1,2,3,4})
   - RK4 is the default (4 stages, 4th order)
   - Tableau-driven implementation

2. **LSERK** (Low-Storage Explicit RK)
   - Reduced memory: only 2 field vectors instead of p+1
   - Useful for large-scale problems
   - Same order as ERK for same stage count

3. **SSPRK** (Strong Stability Preserving RK)
   - TVD property: preserves monotonicity under appropriate CFL restriction
   - SSPRK3(3) is the most common 3-stage, 3rd-order variant

All integrators:
- Advance state `u` from time `t` by step `dt`
- Accept a user-provided RHS callback
- Manage internal stage storage and temporal state


## IO.hpp

VTK input/output utilities for visualizing DG solutions.

**Main functions:**

1. **writeToVTK**
   ```cpp
   void writeToVTK(const Mesh& mesh, const ReferenceElement& ref,
                   const FieldVector& u, const std::string& filename,
                   const std::string& fieldName = "u", Real t = 0)
   ```
   - Writes DG solution to a legacy ASCII VTK file
   - Creates K*Np points (one per DG node) and K cells (one POLY_LINE per element)
   - Paraview visualizes element-local polynomials with discontinuities at interfaces
   - Includes POINT_DATA scalar array with field values

2. **writeTimeSeriesVTK**
   ```cpp
   void writeTimeSeriesVTK(const Mesh& mesh, const ReferenceElement& ref,
                           const FieldVector& u, const std::string& base_name,
                           int frameIndex, Real t, const std::string& fieldName = "u")
   ```
   - Generates time-series VTK files with naming convention: `{base_name}_{frameIndex:06d}.vtk`
   - Useful for animation sequences in Paraview

**Output format:** Legacy ASCII VTK (version 3.0) compatible with Paraview and VisIt.


# Testing

The project includes a comprehensive test suite covering:
- **FieldVector**: storage, indexing, and deep copy
- **ReferenceElement**: GLL node computation, matrix construction, Legendre polynomials
- **Mesh**: uniform and VTK mesh creation, geometry queries, neighbor lookups
- **NumericalFlux**: flux computations for various schemes
- **Integrator**: time stepping with all RK variants
- **Operator**: RHS assembly and interface term handling

Run tests via:
```bash
ctest --test-dir build --output-on-failure
```

Or run individual test binaries directly:
```bash
./build/tests/test_reference_element --gtest_filter="*DMatrix*"
```

# Documentation

Full API documentation is generated by Doxygen. To build:
```bash
doxygen Doxyfile
open docs/html/index.html
```

The documentation describes all classes, methods, and includes algorithm references for key operations (e.g., GLL node computation, Roe flux).

# Compilation notes
**Configure with tests enabled (Debug only)**
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
```
**Build everything**
```bash
cmake --build build --parallel
```

**Run all tests**
```bash
ctest --test-dir build --output-on-failure
```

**Run just one test binary directly (shows individual test names)**
```bash
./build/tests/test_reference_element
```

**Run tests matching a pattern**
```bash
./build/tests/test_reference_element --gtest_filter="*DMatrix*"
```

**Verbose output**
```bash
ctest --test-dir build --output-on-failure -V
```
