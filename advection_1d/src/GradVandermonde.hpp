#pragma once
#include "GradJacobiP.hpp"

template <class RealType, class VectorType, class MatrixType>
MatrixType gradVandermonde(int n, VectorType nodes)
{
  MatrixType grad_vandermonde;
  for (int i = 0; i < nodes.length(); i++)
  {
    for (int j = 0; j < n; j++)
    {
      grad_vandermonde[i][j] = GradJacobiP_point(j,0,0,nodes[i]);
    }
  }
  return grad_vandermonde;
}


