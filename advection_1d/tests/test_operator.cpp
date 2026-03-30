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

    DG::Mesh<double>             mesh  = DG::Mesh<double>::uniform(0.0, 2*M_PI, K);
    DG::ReferenceElement<double> ref   = DG::ReferenceElement<double>(N);
    DG::UpwindFlux<double>       flux  = DG::UpwindFlux<double>(a);
    DG::Operator<double>         op    = DG::Operator<double>(
                                             mesh, ref, flux,
                                             [](double u){ return u; }, a);
    int Np = ref.numDOF();
};

// ── For a constant u=C, rhs must be zero (du_const/dt = 0) ───────────────────
TEST_F(OperatorTest, ConstantSolutionHasZeroRHS)
{
    DG::FieldVector<double> u(K, Np), rhs(K, Np);
    for (int k = 0; k < K; ++k)
        for (int i = 0; i < Np; ++i)
            u.elementPtr(k)[i] = 3.14;

    op.computeRHS(u, rhs);

    for (int k = 0; k < K; ++k)
        for (int i = 0; i < Np; ++i)
            EXPECT_NEAR(rhs.elementPtr(k)[i], 0.0, TOL_TIGHT)
                << "RHS != 0 for constant u at k=" << k << " i=" << i;
}

// ── Energy must not increase at t=0 (dE/dt <= 0) ─────────────────────────────
TEST_F(OperatorTest, EnergyIsNonIncreasing)
{
    DG::FieldVector<double> u(K, Np), rhs(K, Np);

    // Set u = sin(x)
    for (int k = 0; k < K; ++k) {
        double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
        for (int i = 0; i < Np; ++i) {
            double r = ref.nodes()[i];
            u.elementPtr(k)[i] = std::sin(xL + (r+1)*0.5*h);
        }
    }

    op.computeRHS(u, rhs);

    double dEdt = 0;
    for (int k = 0; k < K; ++k) {
        double J = mesh.jacobian(k);
        for (int i = 0; i < Np; ++i)
            dEdt += ref.weights()[i] * J
                  * u.elementPtr(k)[i] * rhs.elementPtr(k)[i];
    }
    dEdt *= 2;
    EXPECT_LE(dEdt, 1e-10) << "dE/dt = " << dEdt << " > 0, scheme is unstable";
}

// ── Continuous solution: flux jumps must be zero at t=0 ──────────────────────
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

// ── RHS for sin(x) should equal -cos(x)/J (exact derivative) ─────────────────
TEST_F(OperatorTest, RHSMatchesExactDerivativeForSmoothSolution)
{
    DG::FieldVector<double> u(K, Np), rhs(K, Np);
    for (int k = 0; k < K; ++k) {
        double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
        for (int i = 0; i < Np; ++i) {
            double r = ref.nodes()[i];
            u.elementPtr(k)[i] = std::sin(xL + (r+1)*0.5*h);
        }
    }

    op.computeRHS(u, rhs);

    // For a continuous solution, rhs = -du/dx = -cos(x)
    // Check only interior nodes (not face nodes which have LIFT correction)
    for (int k = 0; k < K; ++k) {
        double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
        for (int i = 1; i < Np-1; ++i) {   // skip face nodes
            double r     = ref.nodes()[i];
            double x     = xL + (r+1)*0.5*h;
            double exact = -std::cos(x);
            EXPECT_NEAR(rhs.elementPtr(k)[i], exact, TOL_LOOSE)
                << "RHS mismatch at k=" << k << " i=" << i;
        }
    }
}
