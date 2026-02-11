#pragma once

template <class RealType>
RealType JacobiP_point(int n, const RealType a, const RealType b, RealType& x)
{
  RealType p_x_n, a_n, b_n, c_n, d_n;
  RealType p_x_0 = 1.f;
  RealType p_x_1 = (a-b+(a+b+2)*x)/2;
  if (n==0)
  {
    return p_x_0;
  }
  else if (n==1)
  {
    return p_x_1;
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      a_n = 2*(n+1)*(n+a+b+1)*(2*n+a+b);
      b_n = (2*n+a+b+1)*(a*a-b*b);
      c_n = (2*n+a+b)*(2*n+a+b+1)*(2*n+a+b+2);
      d_n = 2*(n+a)*(n+b)*(2*n+a+b+2);

      p_x_n = ((b_n+c_n*x)*P_x_1-d_n*P_x_0)/a_n;

      p_x_0 = p_x_1;
      p_x_1 = p_x_n;
    }
  }
  return p_x_n;
}
