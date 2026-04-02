#pragma once

#include <Eigen/Core>
#include <TNL/Allocators/Default.h>
#include <TNL/Devices/Host.h>
#include <TNL/Matrices/DenseOperations.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Matrices/DenseMatrix.h>
#include <TNL/Matrices/MatrixBase.h>
#include <cmath>
#include <Eigen/Dense>
#include <exception>
#include <iostream>
#include <stdexcept>

namespace DG {

// Gauss-Lobatto nodes and weights (N+1 points, exact for poly deg 2N-1).
template
<
  class Real = double,
  class Index = int
>
class ReferenceElement
{
public:
  using Vector = TNL::Containers::Vector< Real, TNL::Devices::Host, Index >;
  using Matrix = TNL::Matrices::DenseMatrix< Real, TNL::Devices::Host, Index, TNL::Matrices::ElementsOrganization::RowMajorOrder >;

  // Constructor
  explicit ReferenceElement(Index N) : N_(N), Np_(N + 1)
  {
    r_.setSize(Np_);
    w_.setSize(Np_);
    computeGLL_(r_, w_, N_);  // Gauss–Lobatto–Legendre nodes & weights

    V_ = buildVandermonde_(r_); // V_{ij} = P_j(r_i)
    M_ = massMatrix_(V_);
    Vr_ = buildDerivVandermonde_(r_);
    Dr_ = buildDMatrix_(V_, Vr_, M_); // Dr  Vr * inv(V)
    LIFT_ = buildLIFT_(V_); // LIFT = M^{-1} * E, size Np x 2

  }

  // Getters
  Index order() const { return N_; }
  Index numDOF() const { return Np_; }

  const Vector& nodes() const { return r_; }
  const Vector& weights() const { return w_; }
  const Matrix& V() const { return V_; }
  const Matrix Vinv() const { return invertMatrix_(V_); }
  const Matrix& Dr() const { return Dr_; }
  const Matrix& LIFT() const { return LIFT_; }

  // evaluate the Legendre polynomials at x
  // Bonnet's formula see: https://proofwiki.org/wiki/Bonnet%27s_Recursion_Formula
  static Real legendreP(Index n, Real x)
  {
    if ( n < 0 )
    {
      throw std::invalid_argument("The order of the Legendre polynomial must be non-negative ( n >= 0 )");
    }
    Real p0 = 1.0;
    if (n == 0) return p0;
    Real p1 = x;
    if (n == 1) return p1;
    Real norm, p2;
    for (Index k = 1; k < n; ++k)
    {
      p2 = ((2 * k + 1) * x * p1 - k * p0) / (k + 1);
      p0 = p1;
      p1 = p2;
    }
    return p2;
  }

  // normalize the Legendre polynomial
  // normalization n_n = \sqrt{(2n+1)/2}
  static Real legendrePN(Index n, Real x)
  {
    Real norm = std::sqrt(Real(2*n+1) / 2);
    return norm * legendreP(n,x);
  }

  // Evaluate the derivative of the Legendre derivative at x
  // recursion formula, doesn't work for {-1, 1}
  static Real legendrePDeriv(Index n, Real x)
  {
    if (n == 0) return Real(0);
    if (n == 1) return Real(1);

    if (x >= 1 - 1e-16)
    {
      return Real(n) * (Real(n) + 1) / 2.0;
    }
    else if (x <= -1 + 1e-16)
    {
      if (n % 2 == 0)
      {
        return -1.0 * Real(n) * (Real(n) + 1) / 2.0;
      }
      else
      {
        return Real(n) * (Real(n) + 1) / 2.0;
      }
    }

    return n*(legendreP(n-1,x) - x*legendreP(n,x))/(1-x*x);
  }

  // normalize the derivative of the Legendre polynomial
  // normalization n_n = \sqrt{(2n+1)/2}
  static Real legendrePDerivN(Index n, Real x)
  {
    Real norm = std::sqrt(Real(2*n+1) / 2);
    return norm * legendrePDeriv(n,x);
  }

  // Evalueate the second derivative of the Legendre polynomial
  // recursion formula, doesn't work for {-1, 1}
  // add source
  static Real legendrePDeriv2(Index n, Real x)
  {
    if (std::abs(x) >= 1 - 1e-16)
    {
      throw std::invalid_argument("Division by zero in legendrePDeriv2()");
    }

    if (n < 2) return Real(0);
    return (- n * ( n + 1 ) * legendreP(n,x) + 2 * x * legendrePDeriv(n,x))/(1 - x*x);
  }

  // evalueate the third derivative of the Legendre polynomial
  // recursion formula, doesn't work for {-1, 1}
  // add source
  static Real legendrePDeriv3(Index n, Real x)
  {
    if (std::abs(x) >= 1 - 1e-16)
    {
      throw std::invalid_argument("Division by zero in legendrePDeriv3()");
    }

    if (n < 3) return Real(0);
    return (4*x*legendrePDeriv2(n, x)-(n*(n+1)-2)*legendrePDeriv(n,x))/(1-x*x);
  }

  // function for debugging
  void compute_printGLL(Vector& r, Vector& w, const int n)
  {
    for (int i = 2; i < n; i++)
    {
      std::cout << "Nodes and weights for order n = " << i << std::endl;
      r.setSize(i+1);
      w.setSize(i+1);
      computeGLL_(r,w,i);
      std::cout << r << std::endl;
      std::cout << w << std::endl;
      std::cout << std::endl;
    }
  }

private:
  // Gauss–Lobatto–Legendre nodes: endpoints +-1 and interior zeros of dP_N(x)
  // TODO: hardcode the points for low N
  // @param N   polynomial order of the approximation
  void computeGLL_(Vector& r, Vector& w, const int N)
  {
    if (N < 1)
    {
      throw std::invalid_argument("Polynomial order must be at least 1");
    }
    else if ( r.getSize() != N + 1 || w.getSize() != N + 1)
    {
      throw std::invalid_argument("Invalid size of array r || w in function computeGLL_()");
    }

    // endpoints
    r[0] = Real(-1);
    r[N] = Real(+1);
    w[0] = Real(2) / (N * (N + 1));
    w[N] = w[0];


    // TODO: complete this
    // also verify this!
    // switch (n)
    // {
    //   case 3:
    //     r[1] = 0.0;
    //     w[1] = 4.0/3.0;
    //     break;
    // }

    // TODO: verify this code
    // Interior nodes via Newton-Raphson method
    for (Index k = 1; k < N ; k++)
    {
      // Initial guess: Chebyshev-like
      constexpr Real pi = std::acos(Real(-1));
      Real x = - std::cos(pi * k / N);
      for (int iter = 0; iter < 30; iter++)
      {
        Real dPn   = legendrePDeriv(N, x);
        Real ddPn  = legendrePDeriv2(N, x);
        Real dx = dPn / ddPn;
        x -= dx;
        if (std::abs(dx) < 1e-15) break;
      }

      r[k] = x;
      Real Pn = legendreP(N, x);
      w[k] = Real(2) / (N * (N + 1) * Pn * Pn);
    }
  }

  // Vandermonde matrix V_{ij} = sqrt((2j+1)/2) * P_j(r_i)
  // (orthonormal Legendre basis, as in Hesthaven eq. 3.1)
  Matrix buildVandermonde_(const Vector& r) const
  {
    Matrix V(Np_, Np_);
    V.setDimensions(Np_, Np_);
    for (Index i = 0; i < Np_; ++i)
    {
      for (Index j = 0; j < Np_; ++j)
      {
        V(i,j) = legendrePN(j, r[i]);
      }
    }
    return V;
  }

  // derivative of the Vandermonde matrix (normalized)
  // dV_{ij} = sqrt((2j+1)/2) * P'_j(r_i)
  Matrix buildDerivVandermonde_(const Vector& r) const
  {
    Matrix dV(Np_, Np_);
    for (Index i = 0; i < Np_; ++i)
    {
      for (Index j = 0; j < Np_; ++j)
      {
        dV(i,j) = legendrePDerivN(j, r[i]);
      }
    }
    return dV;
  }

  // Convert TNL DenseMatrix to Eigen MatrixXd
  Eigen::MatrixXd tnlToEigen(const Matrix& A) const
  {
    int rows = A.getRows();
    int cols = A.getColumns();
    Eigen::MatrixXd E(rows, cols);
    for (int i = 0; i < rows; i++)
    {
      for (int j = 0; j < cols; j++)
      {
        E(i, j) = A.getElement(i, j);
      }
    }
    return E;
  }

  // Convert Eigen MatrixXd to TNL DenseMatrix
  Matrix eigenToTnl(const Eigen::MatrixXd& E) const
  {
    Matrix A;
    A.setDimensions(E.rows(), E.cols());   // must come before any setElement
    for (int i = 0; i < E.rows(); i++)
    {
      for (int j = 0; j < E.cols(); j++)
      {
        A.setElement(i, j, E(i, j));
      }
    }
    return A;
  }

  // Compute the inverse of a square DenseMatrix using Eigen
  Matrix invertMatrix_(const Matrix& A) const
  {
    Eigen::MatrixXd E = tnlToEigen(A);
    E.inverse();
    return eigenToTnl(E);
  }

  // build the differentiation matrix for the weak formulation
  // D = V dV^T M
  Matrix buildDMatrix_(const Matrix& V, const Matrix& dV, const Matrix& M) const
  {
    Eigen::MatrixXd eV  = tnlToEigen(V);
    Eigen::MatrixXd eM  = tnlToEigen(M);
    Eigen::MatrixXd edV = tnlToEigen(dV);

    Eigen::MatrixXd eD = eV * edV.transpose() * eM;

    return eigenToTnl(eD);
  }

  Matrix massMatrix_(const Matrix& V) const
  {
    Eigen::MatrixXd eV  = tnlToEigen(V);
    Eigen::MatrixXd eVT = eV.transpose();
    Eigen::MatrixXd eD = eV * eVT;

    return eigenToTnl(eD.inverse());
  }

  // build the LIFT matrix
  Matrix buildLIFT_(const Matrix& V) const
  {
    Matrix VVT(Np_, Np_);
    for (Index i = 0; i < Np_; ++i)
        for (Index j = 0; j < Np_; ++j) {
            Real s = Real(0);
            for (Index k = 0; k < Np_; ++k) s += V(i,k) * V(j,k);
            VVT(i,j) = s;
        }
    Matrix L(Np_, 2);
    for (Index i = 0; i < Np_; ++i) {
        L(i,0) = VVT(i, 0);       // left face node (index 0)
        L(i,1) = -VVT(i, Np_-1);   // right face node (index Np-1)
    }
    return L;
  }

  Index  N_, Np_;
  Vector r_, w_;
  Matrix V_, Vr_, Dr_, LIFT_, M_;
};

} // namespace DG
