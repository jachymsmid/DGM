# DGM solver for 1D advection equation

This solver was written along the book Nodal Discontinuous Galerkin Methods by Hesthaven and Warburton.

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

## DGMesh.hpp
This is a mesh class for the DGM. It is templated with `Real`, `Index` and `Device`. It contains
- `numK_` : number of elements in the mesh
- `vertCoods_` : coordinates of the vertices of the elements
- `mesh_` : the underlying TNL mesh

### Constructors
This class has three constructors
- Constructor from an already built mesh
```cpp
DG1DMesh(MeshType mesh);
```
- Constructor that creates a mesh with uniform steps
```cpp
DG1DMesh uniform(Real a, Real b, Index K)
```
- Constructor that reads a VTK file a creates the mesh
```cpp
TODO
```

### Public methods
- `numElements()` : returns the number of elements of the mesh
- `numFaces()` : returns the number of faces of all elements (should be numElements + 1 in 1D)
- `leftVertex(Index k)` : returns the coordinate of the left vertex of k-th element
- `rightVertex(Index k)` : returns the coordinate of the right vertex of k-th element
- `elementSize(Index k)` : returns the size of the k-th element
- `jacobian(Index k)` : returns the Jacobian of the k-th element
- `leftCellOfFace(Index f)` : returns the index of the cell to the left of the f-th face, if the face is the left boundary a -1 is returned
- `rightCellOfFace(Index f)` : returns the index of the cell to the right of the f-th face, if the face is the right boundary a -1 is returned
- `leftNormal()` : returns the left normal (its just -1 in 1D)
- `rightNormal()` : returns the right normal (its just 1 in 1D)
- `isBoundaryFace(Index f)` : checks whether the face with index f is a boundary or not
- `faceCoord(Index f)` : returns the coordinate of the f-th face
- `tnlMesh()` : returns the underlying TNL mesh

### Private methods
- `buildDerivedData_()` :
    - sets the number of elements
    - computes the coordinates of the vertices

## ReferenceElement.hpp
This is a class that holds the data about a reference element. It comprises
- `N_` : the polynomial order of the approximation
- `Np_` : number of degrees of freedom
- `r_` : LGL nodes
- `w_` : LGL weights
- `V_` : the Vandermonde matrix
- `Dr_` : the derivative matrix
- `LIFT_` : the lift

### Public methods
- `order()` : returns the polynomial order of the approximation
- `numDOF()` : returns the number of degrees of freedom
- `nodes()` : returns the coordinates of the LGL nodes on the reference element
- `weights()` : returns the coordinates of the LGL weights on the reference element
- `V()` : returns the Vandermonde matrix
- `Dr()` : returns the derivative matrix
- `LIFT()` : returns the LIFT matrix
- `legendreP(Index n, Real x)` : evaluates the Legendre polynomial of order 'n' at coordinate 'x'
- `legendrePderiv(Index n, Real x)` : evaluates the derivative of the Legendre polynomial of order 'n' at coordinate 'x'

### Private methods
- `computeGLL_(Vector& r, Vector& w)` : computes the GLL nodes and weights
- `buildVandermonde_(const Vector& r)` : builds the Vandermonde matrix from the LGL nodes
- `buildDerivVandermonde(const Vector& r)` : builds the derivative of the Vandermonde matrix
- `invertMatrix_(const Matri& A)` : inverts the matrix A using the Gauss-Jordan elimination
- `buildDMatrix_(const Matrix& V)` : builds the derivative matrix
- `buildLIFT_(const Matrix& V)` : builds the LIFT matrix

## FiledVector.hpp
This is a class that stores
- `data_` : vector of the data
- `K_` : number of elements
- `Np_` : number of degrees of freedom

It has two constructors
- default constructor
- constructor if the 'K' and 'Np' are specified
```cpp
FiledVector(Index K, Index Np);
```

### Public methods
- `numElements()` : returns the number of elements
- `numDOF()` : returns the number of degrees of freedom
- `totalSize()` : returns the total number of data points
- `elementPtr(Index k)` : returns a pointer to the first element of the k-th element
- `data()` : returns the data vector
- `copyFrom(const FiledVector& other)` : deep copy

## NumericalFlux.hpp
This is a general struct that computes the numerical flux from $u^+$, $u^-$ and $\mathbf{n}$
There are two numerical fluxes that inherit from the general struct

### Upwind
```math
begin{align*}
f^* = f(u^-) \text{ for } \mathbf{n} > 0\\
f^* = f(u^+) \text{ for } \mathbf{n} < 0\\
end{align*}
```

### Lax-Friedrichs
```math
f^* = 0.5(f(u^-) + f(u^+)) - 0.5 C (u^+ - u^-)
```
## Operator.hpp
This class computes the RHS for the semi-discrete form of the equation. It comprises
- `mesh_` : the mesh class
- `ref_` : the reference element class
- `flux_` : the numerical flux struct
- `physFlux_` : the physical flux function

It has only one member function `computeRHS(const Field& u, Field& res)` that computes the RHS and saves it into 'res'.

## RK4Integrator.hpp
This class features the 4 step order Runge-Kutta methods for numerical integration. It comprises
- `K_` : number of elements
- `Np_` : number of degrees of freedom
- `callback_` : callback function
- `k1_`,`k2_`,`k3_`,`k4_` : the four derivatives
- `tmp_` : temporary field
- `lastStepCount_` : last step count
- `currentTime_` : current time

## IO.hpp
This is a compilation of functions that work as input output for the solver.

- `writeToVTK(const Mesh mesh, const Ref ref, const string filename, const string fieldname, Real t)` : writes the data to a VTK file
- `writeTimeSeriesVTK(const Mesh mesh, const Ref ref, const string base_name, int frameIndex, Real t, const string fieldname)` : outputs multiple VTK files
- `writePadeVTK(...)` and `writePadeTimeSeriesVTK(...)` : output Padé-Legendre reconstructed fields on a finer equidistant grid

## Padé-Legendre postprocessing

Padé-Legendre reconstruction is available as an **optional postprocessing stage** in `src/main.cpp`.

- Toggle with compile-time constant: `ENABLE_PADE_POSTPROCESS`
- Default degrees: diagonal split `L = N / 2`, `M = N - L`
- Output refinement: `PADE_REFINE_FACTOR` points per original DG node
- Output file series: `output/pade_output_XXXXXX.vtk`

When enabled, the solver keeps the standard DG output unchanged and writes one Padé frame for the **final state** (`t = Tf`).

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
