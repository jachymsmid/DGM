#pragma once

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
    computeGLL_(r_, w_, Np_);  // Gauss–Lobatto–Legendre nodes & weights

    V_ = buildVandermonde_(r_); // V_{ij} = P_j(r_i)
    Vr_ = buildDerivVandermonde_(r_);
    Dr_ = buildDMatrix_(V_, Vr_); // Dr  Vr * inv(V)
    LIFT_ = buildLIFT_(V_); // LIFT = M^{-1} * E, size Np x 2

    std::cout << "Vandermonde matrix:" << std::endl;
    std::cout << V_ << std::endl;
    std::cout << "Derivative of the Vandermonde matrix:" << std::endl;
    std::cout << Vr_ << std::endl;
    std::cout << "LIFT matrix:" << std::endl;
    std::cout << LIFT_ << std::endl;
    std::cout << "Derivation matrix:" << std::endl;
    std::cout << Dr_ << std::endl;
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

  // Evaluate normalized Legendre polynomial P_n at x using three term recurrence
  // Bonnet's formula see: https://proofwiki.org/wiki/Bonnet%27s_Recursion_Formula
  // normalization n_n = \sqrt{(2n+1)/2}
  static Real legendreP(Index n, Real x)
  {
    Real p0 = 1.0/std::sqrt(2);
    if (n == 0) return p0;
    Real p1 = std::sqrt(3.0/2.0)*x;
    if (n == 1) return p1;
    Real norm, p2;
    for (Index k = 1; k < n; ++k)
    {
      norm = std::sqrt(Real(2*k+1) / 2);
      p2 = norm * ((2*k+1)*x*p1 - k*p0) / (k+1);
      p0 = p1;
      p1 = p2;
    }
    return p2;
  }

  // evaluate the unnormalized Legendre polynomials at x
  // Bonnet's formula
  static Real legendrePNN(Index n, Real x)
  {
    Real p0 = 1.0;
    if (n == 0) return p0;
    Real p1 = x;
    if (n == 1) return p1;
    Real norm, p2;
    for (Index k = 1; k < n; ++k)
    {
      p2 = ((2*k+1)*x*p1 - k*p0) / (k+1);
      p0 = p1;
      p1 = p2;
    }
    return p2;
  }

  // Evalueate the derivative of P_n at x
  // Recursion formula
  static Real legendrePDeriv(Index n, Real x)
  {
    if (n == 0) return Real(0);
    if (n == 1) return std::sqrt(3.0/2.0);
    Real norm = std::sqrt(Real(2*n+1) / 2);
    return norm * (( n + 1 ) * legendreP(n-1,x) + x * legendrePDeriv(n-1,x));
  }

  static Real legendrePDerivNN(Index n, Real x)
  {
    if (n == 0) return Real(0);
    if (n == 1) return Real(1);
    return n*(legendrePNN(n-1,x) - x*legendrePNN(n,x))/(1-x*x);
  }

  // Evalueate the second derivative of P_n at x
  // Recursion formula
  static Real legendrePDeriv2(Index n, Real x)
  {
    if (n < 2) return Real(0);
    Real norm = std::sqrt(Real(2*n+1) / 2);
    return norm * (- n * ( n + 1 ) * legendreP(n,x) + 2 * x * legendrePDeriv(n,x))/(1 - x*x);
  }

  static Real legendrePDeriv2NN(Index n, Real x)
  {
    if (n < 2) return Real(0);
    return (2*x*legendrePDerivNN(n, x)-n*(n+1)*legendrePNN(n,x))/(1-x*x);
  }
  
  static Real legendrePDeriv3NN(Index n, Real x)
  {
    if (n < 3) return Real(0);
    return (4*x*legendrePDeriv2NN(n, x)-(n*(n+1)-2)*legendrePDerivNN(n,x))/(1-x*x);
  }
  void compute_printGLL(Vector& r, Vector& w, const int n)
  {
    for (int i = 2; i < n; i++)
    {
      r.setSize(i);
      w.setSize(i);
      computeGLL_(r,w,i);
      std::cout << "Nodes and weights for n = " << i << std::endl;
      std::cout << r << std::endl;
      std::cout << w << std::endl;
      std::cout << std::endl;
    }
  }

private:
  // Gauss–Lobatto–Legendre nodes: endpoints ±1 and interior zeros of dP_{N}(x)
  // TODO: hardcode the points for low N
  void computeGLL_(Vector& r, Vector& w, const int n)
  {
    if (n < 2)
    {
      throw std::runtime_error("Polynomial order must be at least 2");
    }
    // endpoints
    r[0]   = Real(-1);
    r[n-1] = Real(+1);
    w[0]   = Real(2) / (n * (n - 1));
    w[n-1] = w[0];

    if (n == 4)
    {
      r[1] = - std::sqrt(1.0/5.0);
      r[3]
    }

    // Interior nodes via Halley's method iteration on P'_{N}(x) = 0
    for (Index k = 1; k < n - 1; k++)
    {
      // Initial guess: Chebyshev-like
      constexpr Real pi = std::acos(Real(-1));
      Real x = - std::cos(pi * k / N_);
      for (int iter = 0; iter < 30; iter++)
      {
        Real Pn   = legendrePDerivNN(n-1, x);
        Real dPn  = legendrePDeriv2NN(n-1, x);
        Real ddPn = legendrePDeriv3NN(n-1, x);
        Real dx   = 2 * Pn * dPn / (2 * dPn * dPn - Pn * ddPn);
        x += dx;
        if (std::abs(dx) < 1e-15) break;
      }

      r[k] = x;
      Real Pn = legendrePNN(n-1, x);
      w[k] = Real(2) / (n * (n - 1) * Pn * Pn);
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
        V(i,j) = legendreP(j, r[i]);
      }
    }
    return V;
  }

  // Derivative Vandermonde: dV_{ij} = sqrt((2j+1)/2) * P'_j(r_i)
  Matrix buildDerivVandermonde_(const Vector& r) const
  {
    Matrix dV(Np_, Np_);
    dV.setDimensions(Np_, Np_);
    for (Index i = 0; i < Np_; ++i)
    {
      for (Index j = 0; j < Np_; ++j)
      {
        dV(i,j) = legendrePDeriv(j, r[i]);
      }
    }
    return dV;
  }

  // helper functions for matrix inversion - not optimal
  // Convert TNL DenseMatrix to Eigen MatrixXd
  Eigen::MatrixXd tnlToEigen(const Matrix& A) const
  {
    int rows = A.getRows();
    int cols = A.getColumns();
    if (rows == 0 || cols == 0)
      throw std::runtime_error("tnlToEigen: matrix has zero dimensions");
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
  // Compute the inverse of a square DenseMatrix
  Matrix invertMatrix_(const Matrix& A) const
  {
    Eigen::MatrixXd E = tnlToEigen(A);
    E.inverse();
    return eigenToTnl(E);
  }

  // D = dV * inv(V)
  Matrix buildDMatrix_(const Matrix& V, const Matrix& dV) const
  {
    Eigen::MatrixXd eV  = tnlToEigen(V);
    Eigen::MatrixXd edV = tnlToEigen(dV);

    // verify V is not degenerate before inverting
    Real detV = eV.determinant();
    std::cout << "det(V) = " << detV << std::endl;

    Eigen::MatrixXd eVinv = eV.inverse();

    Eigen::MatrixXd eD = edV * eVinv;
    return eigenToTnl(eD);
  }

  Matrix buildLIFT_(const Matrix& V) const
  {
    Matrix L;
    L.setDimensions(Np_, 2);
    for (Index i = 0; i < Np_; ++i) {
        L.setElement(i, 0, Real(0));
        L.setElement(i, 1, Real(0));
    }
    L.setElement(0,      0, Real(1));
    L.setElement(Np_-1,  1, Real(1));
    return L;
  }

  Index  N_, Np_;
  Vector r_, w_;
  Matrix V_, Vr_, Dr_, LIFT_;
};

} // namespace DG
