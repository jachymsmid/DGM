#include <gtest/gtest.h>
#include "ReferenceElement.hpp"

namespace DG = TNL::DGM;

// tolerance for floating-point comparisons
static constexpr double TOL = 1e-12;
static constexpr double TOL2 = 1e-4;

class ReferenceElementTest : public ::testing::TestWithParam<int> {};

// GLL weights must sum to 2
TEST_P(ReferenceElementTest, WeightsSumToTwo)
{
    int N = GetParam();
    DG::ReferenceElement<double> ref(N);
    double sum = 0;
    for (int i = 0; i < ref.numDOF(); ++i)
    {
      sum += ref.weights()[i];
    }
    EXPECT_NEAR(sum, 2.0, TOL);
}

// endpoints must be +1 and -1
TEST_P(ReferenceElementTest, EndpointsArePlusMinusOne)
{
    int N = GetParam();
    DG::ReferenceElement<double> ref(N);
    EXPECT_DOUBLE_EQ(ref.nodes()[0],             -1.0);
    EXPECT_DOUBLE_EQ(ref.nodes()[ref.numDOF()-1], 1.0);
}

// Dmatrix is skew antisymmetric
TEST_P(ReferenceElementTest, DMatrixSkewAntisymmetric)
{
  int N = GetParam();
  DG::ReferenceElement<double> ref(N);
  int Np = ref.numDOF();

  for (int i = 0; i < Np/2; i++)
  {
    for (int j = 0; j < Np/2; j++)
    {
      EXPECT_NEAR(ref.Dr()(i,j), - ref.Dr()(N-i,N-j), TOL);
    }
  }

}

// D * ones = 0  (derivative of a constant is zero)
// row sum should be zero
TEST_P(ReferenceElementTest, DMatrixAnnihilatesConstants)
{
    int N = GetParam();
    DG::ReferenceElement<double> ref(N);
    int Np = ref.numDOF();
    for (int i = 0; i < Np; ++i) {
        double s = 0;
        for (int j = 0; j < Np; ++j)
            s += ref.Dr()(i, j) * 2.0;
        EXPECT_NEAR(s, 0.0, TOL)
            << "D*ones != 0 at row " << i << " for N=" << N;
    }
}

// D * r = 1  (derivative of the identity is one)
TEST_P(ReferenceElementTest, DMatrixDifferentiatesLinear)
{
    int N = GetParam();
    DG::ReferenceElement<double> ref(N);
    int Np = ref.numDOF();
    for (int i = 0; i < Np; ++i) {
        double s = 0;
        for (int j = 0; j < Np; ++j)
            s += ref.Dr()(i, j) * ref.nodes()[j];
        EXPECT_NEAR(s, 1.0, TOL)
            << "D*r != 1 at row " << i << " for N=" << N;
    }
}

// D differentiates polynomials exactly up to degree N
TEST_P(ReferenceElementTest, DMatrixExactForPolynomials)
{
  int N = GetParam() + 1;
  DG::ReferenceElement<double> ref(N);
  int Np = ref.numDOF();

  // Test with f(r) = r^2, f'(r) = 2r
  for (int i = 0; i < Np; ++i)
  {
    double s = 0;
    for (int j = 0; j < Np; ++j)
    {
      s += ref.Dr()(i,j) * ref.nodes()[j] * ref.nodes()[j];
    }
    double exact = 2.0 * ref.nodes()[i];
    EXPECT_NEAR(s, exact, TOL2);
  }
}

// run all tests for N = 1, 2, 3, 4, 5
INSTANTIATE_TEST_SUITE_P(Orders, ReferenceElementTest,
  ::testing::Values(1, 2, 3, 4, 5));

// ── order() and numDOF() getters ──────────────────────────────────────────────
TEST(ReferenceElementGetterTest, OrderAndNumDOF)
{
    for (int N : {1, 3, 5})
    {
        DG::ReferenceElement<double> ref(N);
        EXPECT_EQ(ref.order(),  N);
        EXPECT_EQ(ref.numDOF(), N + 1);
    }
}

// ── Constructor error handling ────────────────────────────────────────────────
TEST(ReferenceElementErrorTest, ConstructorThrowsForN0)
{
    EXPECT_THROW(DG::ReferenceElement<double>(0), std::invalid_argument);
}

TEST(ReferenceElementErrorTest, ConstructorThrowsForNegativeN)
{
    EXPECT_THROW(DG::ReferenceElement<double>(-1), std::invalid_argument);
}

// ── GLL nodes properties ──────────────────────────────────────────────────────
TEST(ReferenceElementNodeTest, NodesAreInRange)
{
    DG::ReferenceElement<double> ref(4);
    int Np = ref.numDOF();
    for (int i = 0; i < Np; ++i)
    {
        EXPECT_GE(ref.nodes()[i], -1.0)
            << "Node " << i << " below -1";
        EXPECT_LE(ref.nodes()[i],  1.0)
            << "Node " << i << " above +1";
    }
}

TEST(ReferenceElementNodeTest, NodesAreStrictlyIncreasing)
{
    DG::ReferenceElement<double> ref(5);
    int Np = ref.numDOF();
    for (int i = 0; i < Np - 1; ++i)
        EXPECT_LT(ref.nodes()[i], ref.nodes()[i + 1])
            << "Nodes not sorted at i=" << i;
}

// ── Legendre polynomial: P_0, P_1, P_2 ───────────────────────────────────────
TEST(LegendrePTest, P0IsConstantOne)
{
    for (double x : {-1.0, -0.5, 0.0, 0.5, 1.0})
        EXPECT_NEAR(DG::ReferenceElement<double>::legendreP(0, x), 1.0, TOL);
}

TEST(LegendrePTest, P1IsIdentity)
{
    for (double x : {-1.0, -0.5, 0.0, 0.5, 1.0})
        EXPECT_NEAR(DG::ReferenceElement<double>::legendreP(1, x), x, TOL);
}

TEST(LegendrePTest, P2KnownValues)
{
    // P_2(x) = (3x^2 - 1) / 2
    auto p2_exact = [](double x) { return (3.0*x*x - 1.0) / 2.0; };
    for (double x : {-1.0, -0.5, 0.0, 0.5, 1.0})
        EXPECT_NEAR(DG::ReferenceElement<double>::legendreP(2, x),
                    p2_exact(x), TOL)
            << "P_2 at x=" << x;
}

TEST(LegendrePTest, NegativeOrderThrows)
{
    EXPECT_THROW(DG::ReferenceElement<double>::legendreP(-1, 0.0),
                 std::invalid_argument);
}

// ── Legendre derivative: P_0', P_1' ──────────────────────────────────────────
TEST(LegendrePDerivTest, P0DerivIsZero)
{
    for (double x : {-0.9, -0.5, 0.0, 0.5, 0.9})
        EXPECT_NEAR(DG::ReferenceElement<double>::legendrePDeriv(0, x), 0.0, TOL);
}

TEST(LegendrePDerivTest, P1DerivIsOne)
{
    for (double x : {-0.9, -0.5, 0.0, 0.5, 0.9})
        EXPECT_NEAR(DG::ReferenceElement<double>::legendrePDeriv(1, x), 1.0, TOL);
}

TEST(LegendrePDerivTest, P2DerivAtBoundary)
{
    // P_2'(1)  = n(n+1)/2 = 3;  P_2'(-1) = -n(n+1)/2 = -3 (n=2 even)
    EXPECT_NEAR(DG::ReferenceElement<double>::legendrePDeriv(2, +1.0),  3.0, TOL);
    EXPECT_NEAR(DG::ReferenceElement<double>::legendrePDeriv(2, -1.0), -3.0, TOL);
}

TEST(LegendrePDerivTest, P2DerivInterior)
{
    // P_2'(x) = 3x
    for (double x : {-0.8, -0.3, 0.0, 0.3, 0.8})
        EXPECT_NEAR(DG::ReferenceElement<double>::legendrePDeriv(2, x), 3.0*x, TOL2);
}

// ── Second and third derivative throw at boundary ────────────────────────────
TEST(LegendrePDeriv2Test, ThrowsAtPlusOne)
{
    EXPECT_THROW(DG::ReferenceElement<double>::legendrePDeriv2(3, +1.0),
                 std::invalid_argument);
}

TEST(LegendrePDeriv2Test, ThrowsAtMinusOne)
{
    EXPECT_THROW(DG::ReferenceElement<double>::legendrePDeriv2(3, -1.0),
                 std::invalid_argument);
}

TEST(LegendrePDeriv2Test, ReturnsZeroForN0)
{
    EXPECT_NEAR(DG::ReferenceElement<double>::legendrePDeriv2(0, 0.0), 0.0, TOL);
}

TEST(LegendrePDeriv2Test, ReturnsZeroForN1)
{
    EXPECT_NEAR(DG::ReferenceElement<double>::legendrePDeriv2(1, 0.0), 0.0, TOL);
}

TEST(LegendrePDeriv3Test, ThrowsAtPlusOne)
{
    EXPECT_THROW(DG::ReferenceElement<double>::legendrePDeriv3(4, +1.0),
                 std::invalid_argument);
}

TEST(LegendrePDeriv3Test, ThrowsAtMinusOne)
{
    EXPECT_THROW(DG::ReferenceElement<double>::legendrePDeriv3(4, -1.0),
                 std::invalid_argument);
}

TEST(LegendrePDeriv3Test, ReturnsZeroForN2)
{
    EXPECT_NEAR(DG::ReferenceElement<double>::legendrePDeriv3(2, 0.0), 0.0, TOL);
}

// ── LIFT matrix dimensions ────────────────────────────────────────────────────
TEST(ReferenceElementLIFTTest, LIFTDimensions)
{
    for (int N : {1, 2, 3, 4})
    {
        DG::ReferenceElement<double> ref(N);
        EXPECT_EQ(ref.LIFT().getRows(),    N + 1) << "LIFT rows wrong for N=" << N;
        EXPECT_EQ(ref.LIFT().getColumns(), 2)      << "LIFT cols wrong for N=" << N;
    }
}

// ── Vandermonde: V * V^{-1} ≈ I ──────────────────────────────────────────────
TEST(ReferenceElementVandermondeTest, VTimesVinvIsIdentity)
{
    DG::ReferenceElement<double> ref(3);
    int Np = ref.numDOF();
    const auto& V    = ref.V();
    const auto  Vinv = ref.Vinv();

    for (int i = 0; i < Np; ++i)
    {
        for (int j = 0; j < Np; ++j)
        {
            double s = 0;
            for (int k = 0; k < Np; ++k)
                s += V(i, k) * Vinv(k, j);
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(s, expected, 1e-10)
                << "V*Vinv != I at (" << i << "," << j << ")";
        }
    }
}
