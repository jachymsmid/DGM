#pragma once
#include "JacobiP.hpp"
#include <cmath>

template <class RealType>
RealType GradJacobiP_point(int n,const RealType alpha, const RealType beta, RealType& x)
{
  RealType b = JacobiP_point(n-1, alpha+1,beta+1,x);
  RealType a = std::sqrt(n*(n+alpha+beta+1))*b

  return a;
}
