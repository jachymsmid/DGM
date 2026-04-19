#include <gtest/gtest.h>
#include "ReferenceElement.hpp"

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
