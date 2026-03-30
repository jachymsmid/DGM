#include <gtest/gtest.h>
#include "../headers/ReferenceElement.hpp"

// Tolerance for floating-point comparisons
static constexpr double TOL = 1e-12;

class ReferenceElementTest : public ::testing::TestWithParam<int> {};

// ── GLL weights must sum to 2 (length of [-1,1]) ─────────────────────────────
TEST_P(ReferenceElementTest, WeightsSumToTwo)
{
    int N = GetParam();
    DG::ReferenceElement<double> ref(N);
    double sum = 0;
    for (int i = 0; i < ref.numDOF(); ++i)
        sum += ref.weights()[i];
    EXPECT_NEAR(sum, 2.0, TOL);
}

// ── Endpoints must be exactly ±1 ─────────────────────────────────────────────
TEST_P(ReferenceElementTest, EndpointsAreMinusOneAndOne)
{
    int N = GetParam();
    DG::ReferenceElement<double> ref(N);
    EXPECT_DOUBLE_EQ(ref.nodes()[0],             -1.0);
    EXPECT_DOUBLE_EQ(ref.nodes()[ref.numDOF()-1], 1.0);
}

// ── D * ones = 0  (derivative of a constant is zero) ─────────────────────────
TEST_P(ReferenceElementTest, DMatrixAnnihilatesConstants)
{
    int N = GetParam();
    DG::ReferenceElement<double> ref(N);
    int Np = ref.numDOF();
    for (int i = 0; i < Np; ++i) {
        double s = 0;
        for (int j = 0; j < Np; ++j)
            s += ref.Dr()(i, j) * 1.0;
        EXPECT_NEAR(s, 0.0, 1e-11)
            << "D*ones != 0 at row " << i << " for N=" << N;
    }
}

// ── D * r = 1  (derivative of the identity is one) ───────────────────────────
TEST_P(ReferenceElementTest, DMatrixDifferentiatesLinear)
{
    int N = GetParam();
    DG::ReferenceElement<double> ref(N);
    int Np = ref.numDOF();
    for (int i = 0; i < Np; ++i) {
        double s = 0;
        for (int j = 0; j < Np; ++j)
            s += ref.Dr()(i, j) * ref.nodes()[j];
        EXPECT_NEAR(s, 1.0, 1e-11)
            << "D*r != 1 at row " << i << " for N=" << N;
    }
}

// ── LIFT must be diag(1/w) at face nodes, zero elsewhere ─────────────────────
TEST_P(ReferenceElementTest, LIFTHasCorrectStructure)
{
    int N = GetParam();
    DG::ReferenceElement<double> ref(N);
    int Np = ref.numDOF();

    // Col 0 (left face): only row 0 is non-zero
    EXPECT_NEAR(ref.LIFT()(0, 0), 1.0 / ref.weights()[0], TOL);
    for (int i = 1; i < Np; ++i)
        EXPECT_NEAR(ref.LIFT()(i, 0), 0.0, TOL)
            << "LIFT(" << i << ",0) should be 0";

    // Col 1 (right face): only row Np-1 is non-zero
    EXPECT_NEAR(ref.LIFT()(Np-1, 1), 1.0 / ref.weights()[Np-1], TOL);
    for (int i = 0; i < Np-1; ++i)
        EXPECT_NEAR(ref.LIFT()(i, 1), 0.0, TOL)
            << "LIFT(" << i << ",1) should be 0";
}

// ── D differentiates polynomials exactly up to degree N ──────────────────────
TEST_P(ReferenceElementTest, DMatrixExactForPolynomials)
{
    int N = GetParam();
    DG::ReferenceElement<double> ref(N);
    int Np = ref.numDOF();

    // Test with f(r) = r^2, f'(r) = 2r
    for (int i = 0; i < Np; ++i) {
        double s = 0;
        for (int j = 0; j < Np; ++j)
            s += ref.Dr()(i,j) * ref.nodes()[j] * ref.nodes()[j];
        double exact = 2.0 * ref.nodes()[i];
        EXPECT_NEAR(s, exact, 1e-10)
            << "D*(r^2) != 2r at node " << i << " for N=" << N;
    }
}

// Run all tests for N = 1, 2, 3, 4, 5
INSTANTIATE_TEST_SUITE_P(Orders, ReferenceElementTest,
    ::testing::Values(1, 2, 3, 4, 5));
