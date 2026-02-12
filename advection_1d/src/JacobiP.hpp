#pragma once

template <class RealType>
RealType JacobiP_point(int N, const RealType a, const RealType b, RealType& x)
{
  RealType p_x_n, a_n, b_n, c_n, d_n, n_r;
  RealType p_x_0 = 1.f;
  RealType p_x_1 = (a-b+(a+b+2)*x)/2;
  if (N==0)
  {
    return p_x_0;
  }
  else if (N==1)
  {
    return p_x_1;
  }
  else
  {
    for (int n = 2; n < N+1; n++)
    {
      n_r = (RealType) n;
      a_n = 2*(n_r+1)*(n_r+a+b+1)*(n_r*2+a+b);
      b_n = (n_r*2+a+b+1)*(a*a-b*b);
      c_n = (2*n_r+a+b)*(2*n_r+a+b+1)*(2*n_r+a+b+2);
      d_n = 2*(n_r+a)*(n_r+b)*(2*n_r+a+b+2);

      p_x_n = ((b_n+c_n*x)*p_x_1-d_n*p_x_0)/a_n;

      p_x_0 = p_x_1;
      p_x_1 = p_x_n;
    }
  }
  return p_x_n;
}
