#pragma once

#include "JacobiGaussQuadrature.hpp"

template <class RealType, class VectorType>
void computeGLnodes(int n, RealType alpha, RealType beta, VectorType& nodes)
{
  QuadratureRule rule = computeGaussJacobi(n-2, alpha+1, beta+1);
  VectorType gq_nodes = rule.nodes;

  for (int i = 1; i < n-1; i++)
  {
    nodes[i] = gq_nodes[i];
  }

  nodes[0] = -1.f;
  nodes[n-1] = 1.f;
}
