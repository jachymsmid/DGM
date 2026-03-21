#pragma once
#include "JacobiGLNodes.hpp"
#include "MatrixDivision.hpp"
#include "DiffMatrix.hpp"

template
<
  class RealType,
  class VectorType,
  class MatrixType,
  int number_of_elements,
  int degree_of_freedom
>
class Mesh
{
public:
  // constructor
  Mesh(MatrixType EToV, VectorType vertex_cord)
  {
    // map the reference elements to physical coordinates
    VectorType va, vb, r;
    for (int i = 0; i < number_of_elements + 1; i++)
    {
      va[i] = EToV[i][0];
      vb[i] = EToV[i][1];
    }

    // LGL nodes on a reference element
    r.reserve(degree_of_freedom);
    computeGLnodes(degree_of_freedom, 0.f, 0.f , r); 

    // physical coordinates of all the nodes
    for (int i = 0; i < number_of_elements; i++)
    {
      for (int j = 0; j < degree_of_freedom; j++)
      {
        _x_cord[i][j] = vertex_cord[va[i]] + 0.5 * (r[j] + 1) * (vertex_cord[vb[i]] - vertex_cord[va[i]]);
      }
    }

    // compute the metric of the mapping
    for (int i = 0; i < degree_of_freedom; i++)
    {
      // cant really figure out the dimensions of the result?
    }

    // physical coordinates of the vertices
    for (int i = 0; i < number_of_elements; i++)
    {
      // is this really needed?
    }

    // normals for each element
    for (int i = 0; i < number_of_elements; i++)
    {
      _normals[i][0] = -1.f;
      _normals[i][1] = 1.f;
    }

    // connection
    _FToV = ;

    // face connections
    _FToF = matrixMultiplication(FToV, FToV.transpose());
    // remove self reference
    for (int i = 0; i < 2*number_of_elements; i++)
    {
      _FToF[i][i] -= 1;
    }

  }

private:
  VectorType _vertex_cord;
  MatrixType _EToV;
  VectorType _x_cord;
  int _number_elements;
  int _dof;
  VectorType _rx;
  VectorType _jacobian;
  MatrixType _normals;
  MatrixType _FToF;
  MatrixType _FToV;
};
