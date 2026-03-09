#pragma once

template <class RealType, class VectorType, class FunctionType>
void rk4(VectorType& data, RealType t, RealType delta_t, FunctionType rhs)
{

  for (int i = 0; i < data.size(); i++)
  {
    RealType k1 = rhs(data[i], t);
    RealType k2 = rhs(data[i] + 0.5*delta_t*k1, t + delta_t);
    RealType k3 = rhs(data[i] + 0.5*delta_t*k2, t + delta_t);
    RealType k4 = rhs(data[i] + 0.5*delta_t*k3, t + delta_t);
    data[i] = 1/6*delta_t*(k1[i]+k2[i]+k3[i]+k4[i]);
  }
}

