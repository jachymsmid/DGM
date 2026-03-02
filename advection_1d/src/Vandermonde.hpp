#pragma once
#include "JacobiP.hpp"

template <class RealType, class VectorType, class MatrixType>
void computeVandermonde(int n, VectorType nodes)
{
  MatrixType vandermonde;
  for (int i = 0; i < nodes.length(); i++)
  {
    for (int j = 0; j < n; j++)
    {
      vandermonde[i][j] = JacobiP_point(nodes[i], 0, 0, j);
    }
  }
  return vandermonde;
}
