#pragma once

#include <TNL/Allocators/Default.h>
#include <TNL/Devices/Host.h>
#include <TNL/Matrices/DenseOperations.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Matrices/DenseMatrix.h>
#include <TNL/Matrices/MatrixBase.h>
#include <cmath>
#include <Eigen/Dense>

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
    computeGLL_(r_, w_);  // Gauss–Lobatto–Legendre nodes & weights

    V_ = buildVandermonde_(r_); // V_{ij} = P_j(r_i)
    Dr_ = buildDMatrix_(V_); // Dr = V * (dV/dr)^{-1} ... actually D = Vr * inv(V)
    LIFT_ = buildLIFT_(V_); // LIFT = M^{-1} * E, size Np x 2
  }

  // Getters
  Index order() const { return N_; }
  Index numDOF() const { return Np_; }

  const Vector& nodes() const { return r_; }
  const Vector& weights() const { return w_; }
  const Matrix& V() const { return V_; }
  const Matrix& Dr() const { return Dr_; }
  const Matrix& LIFT() const { return LIFT_; }

  // Evaluate Legendre polynomial P_n at x using three term recurrence
  // Bonnet's formula see: https://proofwiki.org/wiki/Bonnet%27s_Recursion_Formula
  static Real legendreP(Index n, Real x)
  {
    if (n == 0) return Real(1);
    if (n == 1) return x;
    Real p0 = 1, p1 = x, p2 = 0;
    for (Index k = 1; k < n; ++k)
    {
      p2 = ((2*k+1)*x*p1 - k*p0) / (k+1);
      p0 = p1; p1 = p2;
    }
    return p2;
  }

  // Evalueate the derivative of P_n at x
  // Recursion formula
  static Real legendrePDeriv(Index n, Real x)
  {
    if (n == 0) return Real(0);
    if (n == 1) return Real(1);
    return ( n + 1 ) * legendreP(n,x) + x * legendrePDeriv(n,x);
  }

private:
  // Gauss–Lobatto–Legendre nodes: endpoints ±1 and interior zeros of dP_{N}(x)
  // TODO: hardcode the points for low N
  void computeGLL_(Vector& r, Vector& w)
  {
    const int n = Np_;
    // endpoints
    r[0]   = Real(-1);
    r[n-1] = Real(+1);
    w[0]   = Real(2) / (N_ * (N_ + 1));
    w[n-1] = w[0];

    // Interior nodes via Newton iteration on P'_{N}(x) = 0
    for (Index k = 1; k < n - 1; k++)
    {
      // Initial guess: Chebyshev-like
      constexpr Real pi = std::acos(Real(-1));
      Real x = - std::cos(pi * k / N_);
      for (int iter = 0; iter < 30; iter++)
      {
        // P'_N and P''_N
        Real Pn   = legendreP(N_, x);
        Real Pnm1 = legendreP(N_-1, x);
        // P'_N(x) = N*(P_{N-1}(x) - x*P_N(x))/(1-x^2)
        Real dPn  = N_ * (Pnm1 - x * Pn) / (1 - x*x + 1e-16);
        Real ddPn = (2*x*dPn - N_*(N_+1)*Pn) / (1 - x*x + 1e-16);
        Real dx   = -dPn / ddPn;
        x += dx;
        if (std::abs(dx) < 1e-15) break;
      }

      r[k] = x;
      Real Pn = legendreP(N_, x);
      w[k] = Real(2) / (N_ * (N_+1) * Pn * Pn);
    }
  }

  // Vandermonde matrix V_{ij} = sqrt((2j+1)/2) * P_j(r_i)
  // (orthonormal Legendre basis, as in Hesthaven eq. 3.1)
  Matrix buildVandermonde_(const Vector& r) const
  {
    Matrix V(Np_, Np_);
    for (Index i = 0; i < Np_; ++i)
    {
      for (Index j = 0; j < Np_; ++j)
      {
        Real norm = std::sqrt(Real(2*j+1) / 2);
        V(i,j) = norm * legendreP(j, r[i]);
      }
    }
    return V;
  }

  // Derivative Vandermonde: dV_{ij} = sqrt((2j+1)/2) * P'_j(r_i)
  Matrix buildDerivVandermonde_(const Vector& r) const
  {
    Matrix dV(Np_, Np_);
    for (Index i = 0; i < Np_; ++i)
    {
      for (Index j = 0; j < Np_; ++j)
      {
        Real norm = std::sqrt(Real(2*j+1) / 2);
        dV(i,j) = norm * legendrePDeriv(j, r[i]);
      }
    }
    return dV;
  }

  // helper functions for matrix inversion - not optimal
  // Convert TNL DenseMatrix to Eigen MatrixXd
  Eigen::MatrixXd tnlToEigen(const Matrix& A)
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
  Matrix eigenToTnl(const Eigen::MatrixXd& E)
  {
    Matrix A(E.rows(), E.cols());
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
  Matrix buildDMatrix_(const Matrix& V) const
  {
    Matrix dV  = buildDerivVandermonde_(r_);
    Matrix Vinv = invertMatrix_(V);

    // D = dV * Vinv
    Matrix D(Np_, Np_);
    for (Index i = 0; i < Np_; ++i)
    {
      for (Index j = 0; j < Np_; ++j)
      {
        Real s = Real(0);
        for (Index k = 0; k < Np_; ++k)
        {
          s += dV(i, k) * Vinv(k, j);
        }
        D.setElement(i, j, s);
      }
    }
    return D;
  }

  Matrix buildLIFT_(const Matrix& V) const
  {
    Matrix E(Np_, 2);
    Matrix L(Np_, Np_);
    for (Index i = 1; i < Np_ - 1; i++)
    {
      E.setElement(i, 0, Real(0));
      E.setElement(i, 1, Real(0));
    }

    E.setElement(0, 0, Real(1));
    E.setElement(Np_ - 1, Np_ - 1, Real(1));
    getMatrixProduct(L, V.getInPlaceTransposition(), E);
    getMatrixProduct(L , V, L);
    return L;
  }

  Index  N_, Np_;
  Vector r_, w_;
  Matrix V_, Dr_, LIFT_;
};

} // namespace DG
