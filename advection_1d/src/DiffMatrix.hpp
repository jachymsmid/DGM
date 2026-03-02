#pragma once
#include "MatrixDivision.hpp"
#include "GradVandermonde.hpp"

template <class VectorType, class MatrixType>
MatrixType DiffMatrix(int n, VectorType nodes, MatrixType V)
{
  MatrixType b = gradVandermonde(n, nodes);

  MatrixType diffMatrix = divideMatrices(b,V);

  return diffMatrix;
}
