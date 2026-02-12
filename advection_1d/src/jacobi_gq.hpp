
void jacobi_compute ( int order, double alpha, double beta, double xtab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_COMPUTE computes a Gauss-Jacobi quadrature rule.
//
//  Discussion:
//
//    The integration interval is [ -1, 1 ].
//
//    The weight function is w(x) = (1-X)^ALPHA * (1+X)^BETA.
//
//    The integral to approximate:
//
//      Integral ( -1 <= X <= 1 ) (1-X)^ALPHA * (1+X)^BETA * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//    Thanks to Xu Xiang of Fudan University for pointing out that
//    an earlier implementation of this routine was incorrect!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest
//    C++ version by John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the quadrature rule to be computed.
//
//    Input, double ALPHA, BETA, the exponents of (1-X) and
//    (1+X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
//    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  double an;
  double *b;
  double bn;
  double *c;
  double cc;
  double delta;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double r3;
  double temp;
  double x;

  b = new double[order];
  c = new double[order];
//
//  Check ALPHA and BETA.
//
  if ( alpha <= -1.0 )
  {
    cout << "\n";
    cout << "JACOBI_COMPUTE - Fatal error!\n";
    cout << "  -1.0 < ALPHA is required.\n";
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    cout << "\n";
    cout << "JACOBI_COMPUTE - Fatal error!\n";
    cout << "  -1.0 < BETA is required.\n";
    exit ( 1 );
  }
//
//  Set the recursion coefficients.
//
  for ( i = 1; i <= order; i++ )
  {
    if ( alpha + beta == 0.0 || beta - alpha == 0.0 )
    {
      b[i-1] = 0.0;
    }
    else
    {
      b[i-1] = ( alpha + beta ) * ( beta - alpha ) / 
             ( ( alpha + beta + ( double ) ( 2 * i ) ) 
             * ( alpha + beta + ( double ) ( 2 * i - 2 ) ) );
    }

    if ( i == 1 )
    {
      c[i-1] = 0.0;
    }
    else
    {
      c[i-1] = 4.0 * ( double ) ( i - 1 ) 
         * ( alpha + ( double ) ( i - 1 ) ) 
          * ( beta + ( double ) ( i - 1 ) ) 
            * ( alpha + beta + ( double ) ( i - 1 ) ) / 
            ( ( alpha + beta + ( double ) ( 2 * i - 1 ) ) 
            * pow ( alpha + beta + ( double ) ( 2 * i - 2 ), 2 ) 
            * ( alpha + beta + ( double ) ( 2 * i - 3 ) ) );
    }
  }

  delta = r8_gamma ( alpha        + 1.0 ) 
        * r8_gamma (         beta + 1.0 ) 
        / r8_gamma ( alpha + beta + 2.0 );

  prod = 1.0;
  for ( i = 2; i <= order; i++ )
  {
    prod = prod * c[i-1];
  }
  cc = delta * pow ( 2.0, alpha + beta + 1.0 ) * prod;

  for ( i = 1; i <= order; i++ )
  {
    if ( i == 1 )
    {
      an = alpha / ( double ) ( order );
      bn = beta / ( double ) ( order );

      r1 = ( 1.0 + alpha ) 
        * ( 2.78 / ( 4.0 + ( double ) ( order * order ) ) 
        + 0.768 * an / ( double ) ( order ) );

      r2 = 1.0 + 1.48 * an + 0.96 * bn 
        + 0.452 * an * an + 0.83 * an * bn;

      x = ( r2 - r1 ) / r2;
    }
    else if ( i == 2 )
    {
      r1 = ( 4.1 + alpha ) / 
        ( ( 1.0 + alpha ) * ( 1.0 + 0.156 * alpha ) );

      r2 = 1.0 + 0.06 * ( ( double ) ( order ) - 8.0 ) * 
        ( 1.0 + 0.12 * alpha ) / ( double ) ( order );

      r3 = 1.0 + 0.012 * beta * 
        ( 1.0 + 0.25 * r8_abs ( alpha ) ) / ( double ) ( order );

      x = x - r1 * r2 * r3 * ( 1.0 - x );
    }
    else if ( i == 3 )
    {
      r1 = ( 1.67 + 0.28 * alpha ) / ( 1.0 + 0.37 * alpha );

      r2 = 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order );

      r3 = 1.0 + 8.0 * beta / 
        ( ( 6.28 + beta ) * ( double ) ( order * order ) );

      x = x - r1 * r2 * r3 * ( xtab[0] - x );
    }
    else if ( i < order - 1 )
    {
      x = 3.0 * xtab[i-2] - 3.0 * xtab[i-3] + xtab[i-4];
    }
    else if ( i == order - 1 )
    {
      r1 = ( 1.0 + 0.235 * beta ) / ( 0.766 + 0.119 * beta );

      r2 = 1.0 / ( 1.0 + 0.639 
        * ( ( double ) ( order ) - 4.0 ) 
        / ( 1.0 + 0.71 * ( ( double ) ( order ) - 4.0 ) ) );

      r3 = 1.0 / ( 1.0 + 20.0 * alpha / ( ( 7.5 + alpha ) * 
        ( double ) ( order * order ) ) );

      x = x + r1 * r2 * r3 * ( x - xtab[i-3] );
    }
    else if ( i == order )
    {
      r1 = ( 1.0 + 0.37 * beta ) / ( 1.67 + 0.28 * beta );

      r2 = 1.0 / 
        ( 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order ) );

      r3 = 1.0 / ( 1.0 + 8.0 * alpha / 
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) ) );

      x = x + r1 * r2 * r3 * ( x - xtab[i-3] );
    }

    jacobi_root ( &x, order, alpha, beta, &dp2, &p1, b, c );

    xtab[i-1] = x;
    weight[i-1] = cc / ( dp2 * p1 );
  }
//
//  Reverse the order of the values.
//
  for ( i = 1; i <= order/2; i++ )
  {
    temp          = xtab[i-1];
    xtab[i-1]     = xtab[order-i];
    xtab[order-i] = temp;
  }

  for ( i = 1; i <=order/2; i++ )
  {
    temp            = weight[i-1];
    weight[i-1]     = weight[order-i];
    weight[order-i] = temp;
  }

  delete [] b;
  delete [] c;

  return;
}
