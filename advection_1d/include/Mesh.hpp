/**
 * @file Mesh.hpp
 * @brief DG mesh wrapper built on top of TNL meshes for 1D problems.
 *
 * The Mesh class provides element counts, vertex coordinates, element
 * sizes, jacobians and neighbour/face queries tailored for 1D DG
 * discretizations. Factory helpers include uniform(...) and readVTK(...).
 */

#pragma once

#include "MeshConfig.hpp"
#include <TNL/Containers/Array.h>
#include <TNL/Containers/StaticVector.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Meshes/Readers/VTKReader.h>
#include <stdexcept>

namespace TNL::DGM {

/**
 * @class Mesh
 * @brief 1D DG mesh wrapper providing element geometry and adjacency queries.
 *
 * Built on top of TNL meshes; provides convenient accessors for element
 * vertices, sizes, jacobians and face neighbour queries used by the DG
 * operator and I/O routines.
 */
template
<
  class Real = double,
  class Device = TNL::Devices::Host,
  class Index = int
>
class Mesh
{
public:
  using TNLMeshType = TNLMesh<Real, Device, Index>; // my alias from MeshConfig.hpp
  using RealArray = TNL::Containers::Array<Real,  Device, Index>;
  using IndexArray = TNL::Containers::Array<Index, Device, Index>;

  // Sentinel for boundary faces (no right/left neighbour)
  static constexpr Index BOUNDARY_FACE = -1;

  //  Construction from a TNL mesh
  explicit Mesh(TNLMeshType mesh)
      : mesh_(std::move(mesh))
  {
    buildDerivedData_();
  }

  //  Factory: uniform mesh on [a,b] with K elements
  /**
   * @brief Build a uniform mesh on interval [a,b] with K elements.
   * @param a left boundary
   * @param b right boundary
   * @param K number of elements
   * @return Mesh instance
   */
  static Mesh uniform(Real a, Real b, Index K)
  {
    if (K < 1) throw std::invalid_argument("K must be >= 1");
    using Builder = TNL::Meshes::MeshBuilder<TNLMeshType>;
    Builder builder;
    builder.setEntitiesCount(K + 1, K);   // K+1 vertices, K cells

    const Real h = (b - a) / K;

    for (Index i = 0; i <= K; ++i)
      builder.setPoint(i, TNL::Containers::StaticVector<1,Real>{a + i * h});

    for (Index k = 0; k < K; ++k) {
      builder.getCellSeed(k).setCornerId(0, k);
      builder.getCellSeed(k).setCornerId(1, k + 1);
    }
    TNLMeshType m;
    builder.build(m);
    return Mesh(std::move(m));
  }

  // ------------------------------------------------------------------ //
  //  VTK Reader: construct mesh from a VTK file
  // ------------------------------------------------------------------ //

  /**
   * @brief Construct a Mesh by reading a VTK file produced by this code.
   * @param filename path to the VTK mesh file
   * @return Mesh constructed from the file
   * @throws std::runtime_error on unsupported mesh dimension
   */
  static Mesh readVTK(const std::string& filename)
  {
    using TNLMeshType = typename Mesh<Real, Device, Index>::TNLMeshType;

    TNL::Meshes::Readers::VTKReader reader(filename);
    reader.detectMesh();

    // Sanity checks: must be a 1D mesh of line segments
    if (reader.getMeshDimension() != 1)
      throw std::runtime_error(
        "readMeshFromVTK: file '" + filename +
        "' has mesh dimension " +
        std::to_string(reader.getMeshDimension()) +
        ", expected 1.");

    TNLMeshType tnlMesh;
    reader.loadMesh(tnlMesh);

    return Mesh<Real, Device, Index>(std::move(tnlMesh));
  }


  // ------------------------------------------------------------------ //
  //  Getters
  // ------------------------------------------------------------------ //
  /**
   * @brief Number of elements in the mesh (K).
   * @return K
   */
  Index numElements() const { return numK_; }
  Index numFaces()  const { return numK_ + 1; } // K+1 vertices in 1D

  // Physical vertex coordinates of element k: [x_L, x_R]
  /**
   * @brief Coordinate of the left vertex of element k.
   * @param k element index
   * @return x coordinate
   */
  Real leftVertex (Index k) const { return vertCoords_[k];     }
  /**
   * @brief Coordinate of the right vertex of element k.
   * @param k element index
   * @return x coordinate
   */
  Real rightVertex(Index k) const { return vertCoords_[k + 1]; }
  /**
   * @brief Physical length of element k (x_R - x_L).
   * @param k element index
   */
  Real elementSize(Index k) const { return rightVertex(k) - leftVertex(k); }
  Real minElementSize() const { return min_h_; }
  Real maxElementSize() const { return max_h_; }

  // Jacobian of the affine map r -> x:  x = x_L + (r+1)/2 * h_k
  Real jacobian(Index k) const { return elementSize(k) * Real(0.5); }
  // min(J_i)
  Real minJacobian() const { return Real(0.5) * minElementSize(); }
  // max(J_i)
  Real maxJacobian() const { return Real(0.5) * maxElementSize(); }

  // Left/right cell indices sharing face f (BOUNDARY_FACE if none)
  // Face f = vertex f.  Faces 0..K:
  //   face 0 is the left  boundary (only right cell k=0 exists)
  //   face K is the right boundary (only left  cell k=K-1 exists)
  //   face f (0 < f < K): left cell = f-1, right cell = f
  Index leftCellOfFace(Index f) const { return (f == 0) ? BOUNDARY_FACE : f - 1; }
  Index rightCellOfFace(Index f) const { return (f == numK_) ? BOUNDARY_FACE : f; }

  // Outward normal of cell k at its left face (-1) and right face (+1)
  static constexpr Real leftNormal() { return Real(-1); }
  static constexpr Real rightNormal() { return Real(+1); }

  // Is face f on the domain boundary?
  bool isBoundaryFace(Index f) const { return f == 0 || f == numK_; }

  // Physical coordinate of face f
  Real faceCoord(Index f) const { return vertCoords_[f]; }

  // The underlying TNL mesh (for VTK I/O etc.)
  const TNLMeshType& tnlMesh() const { return mesh_; }

private:
  void buildDerivedData_()
  {
    // number of cells in a TNL mesh
    numK_ = mesh_.template getEntitiesCount<1>();

    // number of vertices
    vertCoords_.setSize(numK_ + 1);

    Real h = 0;
    Real vertCoordPrev;

    // for every mesh cell
    // set coordinates of the vertices, get max and min spatial step
    for (Index i = 0; i <= numK_; ++i)
    {
      auto pt = mesh_.template getEntity<0>(i).getPoint();
      vertCoords_[i] = pt[0];
      if (i == 0)
      {
        vertCoordPrev = pt[0];
        continue;
      }
      h = vertCoords_[i] - vertCoordPrev;
      max_h_ = std::max(h, max_h_);
      min_h_ = std::min(h, min_h_);
      vertCoordPrev = pt[0];
    }
  }

  TNLMeshType  mesh_;
  Index numK_{0};
  Real max_h_{0.0};
  Real min_h_{100.0};
  RealArray vertCoords_;   // x_0, x_1, ..., x_K
};

} // namespace DG
