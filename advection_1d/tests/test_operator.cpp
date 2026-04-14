#include <functional>
#include <gtest/gtest.h>
#include "../headers/Mesh.hpp"
#include "../headers/ReferenceElement.hpp"
#include "../headers/FieldVector.hpp"
#include "../headers/NumericalFlux.hpp"
#include "../headers/Operator.hpp"
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
