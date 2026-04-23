#include <functional>
#include <gtest/gtest.h>
#include "Mesh.hpp"
#include "ReferenceElement.hpp"
#include "FieldVector.hpp"
#include "NumericalFlux.hpp"
#include "Operator.hpp"
#include <cmath>

static constexpr double TOL_TIGHT = 1e-11;
static constexpr double TOL_LOOSE = 1e-6;

class OperatorTest : public ::testing::Test
{
protected:
  static constexpr int K = 10, N = 4;
  static constexpr double a = 1.0;

  using Real = double;
  
  std::function<Real(Real)> physical_flux = [&] ( Real u ) -> Real { return a * u; };
  std::function<Real(Real)> advection_speed = [&] ( Real u ) -> Real { return a; };

  DG::Mesh<double> mesh  = DG::Mesh<double>::uniform(0.0, 2*M_PI, K);
  DG::ReferenceElement<double> ref = DG::ReferenceElement<double>(N);
  DG::UpwindFlux<double> flux = DG::UpwindFlux<double>( advection_speed, physical_flux );
  DG::Operator<double> op = DG::Operator<double>(mesh, ref, flux, [](double u){ return u; });
  int Np = ref.numDOF();
};

// for a continuous solution: flux jumps must be zero at t=0
TEST_F(OperatorTest, ContinuousSolutionHasNoFluxJumps)
{
    DG::FieldVector<double> u(K, Np);
    for (int k = 0; k < K; ++k) {
        double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
        for (int i = 0; i < Np; ++i) {
            double r = ref.nodes()[i];
            u.elementPtr(k)[i] = std::sin(xL + (r+1)*0.5*h);
        }
    }

    // All element boundaries should match
    for (int k = 0; k < K-1; ++k) {
        double u_right = u.elementPtr(k)[Np-1];
        double u_left  = u.elementPtr(k+1)[0];
        EXPECT_NEAR(u_right, u_left, TOL_TIGHT)
            << "Discontinuity at face " << k+1;
    }
}

// RHS for sin(x) should equal -cos(x)/J (exact derivative)
// is this even true?
TEST_F(OperatorTest, RHSMatchesExactDerivativeForSmoothSolution)
{
  DG::FieldVector<double> u(K, Np), rhs(K, Np);
  for (int k = 0; k < K; ++k)
  {
    double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
    for (int i = 0; i < Np; ++i)
    {
      double r = ref.nodes()[i];
      u.elementPtr(k)[i] = std::sin(xL + (r+1)*0.5*h);
    }
  }

  op.computeRHS(u, rhs, 1.0);

  // For a continuous solution, rhs = -du/dx = -cos(x)
  // Check only interior nodes (not face nodes which have LIFT correction)
  for (int k = 0; k < K; ++k)
  {
    double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
    for (int i = 1; i < Np-1; ++i) // skip face nodes
    {
      double r     = ref.nodes()[i];
      double x     = xL + (r+1)*0.5*h;
      double exact = -std::cos(x);
      EXPECT_NEAR(rhs.elementPtr(k)[i], exact, TOL_LOOSE)
          << "RHS mismatch at k=" << k << " i=" << i;
    }
  }
}

// ── Zero field yields zero RHS ────────────────────────────────────────────────
// For f(u)=u and u=0, all fluxes are zero, so rhs must be identically zero.
TEST_F(OperatorTest, ZeroFieldYieldsZeroRHS)
{
    DG::FieldVector<double> u(K, Np);   // zero-initialised
    DG::FieldVector<double> rhs(K, Np);
    op.computeRHS(u, rhs, 0.0);

    for (int k = 0; k < K; ++k)
        for (int i = 0; i < Np; ++i)
            EXPECT_NEAR(rhs.elementPtr(k)[i], 0.0, TOL_TIGHT)
                << "RHS non-zero at k=" << k << " i=" << i;
}

// ── Operator produces finite values for LaxFriedrichs flux ───────────────────
TEST_F(OperatorTest, LaxFriedrichsFluxProducesFiniteRHS)
{
    DG::LaxFriedrichsFlux<double> lfFlux(advection_speed, physical_flux);
    DG::Operator<double> lfOp(mesh, ref, lfFlux, physical_flux);

    DG::FieldVector<double> u(K, Np), rhs(K, Np);
    for (int k = 0; k < K; ++k)
    {
        double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
        for (int i = 0; i < Np; ++i)
        {
            double r = ref.nodes()[i];
            u.elementPtr(k)[i] = std::sin(xL + (r + 1.0) * 0.5 * h);
        }
    }

    ASSERT_NO_THROW(lfOp.computeRHS(u, rhs, 0.0));

    for (int k = 0; k < K; ++k)
        for (int i = 0; i < Np; ++i)
            EXPECT_TRUE(std::isfinite(rhs.elementPtr(k)[i]))
                << "Non-finite RHS with LF flux at k=" << k << " i=" << i;
}

// ── Operator produces finite values for Roe flux ─────────────────────────────
TEST_F(OperatorTest, RoeFluxProducesFiniteRHS)
{
    DG::RoeFlux<double> roeFlux(advection_speed, physical_flux);
    DG::Operator<double> roeOp(mesh, ref, roeFlux, physical_flux);

    DG::FieldVector<double> u(K, Np), rhs(K, Np);
    for (int k = 0; k < K; ++k)
    {
        double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
        for (int i = 0; i < Np; ++i)
        {
            double r = ref.nodes()[i];
            u.elementPtr(k)[i] = std::sin(xL + (r + 1.0) * 0.5 * h);
        }
    }

    ASSERT_NO_THROW(roeOp.computeRHS(u, rhs, 0.0));

    for (int k = 0; k < K; ++k)
        for (int i = 0; i < Np; ++i)
            EXPECT_TRUE(std::isfinite(rhs.elementPtr(k)[i]))
                << "Non-finite RHS with Roe flux at k=" << k << " i=" << i;
}

// ── Operator produces finite values for Godunov flux ─────────────────────────
TEST_F(OperatorTest, GodunovFluxProducesFiniteRHS)
{
    DG::GodunovFlux<double> godFlux(advection_speed, physical_flux);
    DG::Operator<double> godOp(mesh, ref, godFlux, physical_flux);

    DG::FieldVector<double> u(K, Np), rhs(K, Np);
    for (int k = 0; k < K; ++k)
    {
        double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
        for (int i = 0; i < Np; ++i)
        {
            double r = ref.nodes()[i];
            u.elementPtr(k)[i] = std::sin(xL + (r + 1.0) * 0.5 * h);
        }
    }

    ASSERT_NO_THROW(godOp.computeRHS(u, rhs, 0.0));

    for (int k = 0; k < K; ++k)
        for (int i = 0; i < Np; ++i)
            EXPECT_TRUE(std::isfinite(rhs.elementPtr(k)[i]))
                << "Non-finite RHS with Godunov flux at k=" << k << " i=" << i;
}

// ── rhsFunction() wrapper produces identical output to computeRHS ────────────
TEST_F(OperatorTest, RhsFunctionMatchesComputeRHS)
{
    DG::FieldVector<double> u(K, Np), rhs1(K, Np), rhs2(K, Np);
    for (int k = 0; k < K; ++k)
    {
        double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
        for (int i = 0; i < Np; ++i)
        {
            double r = ref.nodes()[i];
            u.elementPtr(k)[i] = std::sin(xL + (r + 1.0) * 0.5 * h);
        }
    }

    op.computeRHS(u, rhs1, 0.5);
    auto f = op.rhsFunction();
    f(u, rhs2, 0.5);

    for (int k = 0; k < K; ++k)
        for (int i = 0; i < Np; ++i)
            EXPECT_NEAR(rhs1.elementPtr(k)[i], rhs2.elementPtr(k)[i], TOL_TIGHT)
                << "rhsFunction mismatch at k=" << k << " i=" << i;
}
