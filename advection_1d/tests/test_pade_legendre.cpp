#include <gtest/gtest.h>
#include "IO.hpp"
#include "Mesh.hpp"
#include "Postprocessing.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

namespace DG = TNL::DGM;

// ── Tolerances ────────────────────────────────────────────────────────────────
static constexpr double TOL       = 1e-10;  // tight (analytical)
static constexpr double TOL_LOOSE = 1e-6;   // loose (numerical quadrature)
static constexpr double TOL_PADE  = 1e-4;   // Padé coefficient agreement

// ─────────────────────────────────────────────────────────────────────────────
//  modalCoeffs: constant polynomial u_h = A
// ─────────────────────────────────────────────────────────────────────────────
TEST(ModalCoeffsTest, ConstantPolynomial)
{
    DG::ReferenceElement<double> ref(4);
    DG::PadeLegendreSolver<double> solver(ref, 2, 2);

    const int Np = ref.numDOF();
    std::vector<double> nodal(Np, 3.7);

    const auto c = solver.modalCoeffs(nodal.data());

    EXPECT_NEAR(c[0], 3.7, TOL);
    for (int n = 1; n < Np; ++n)
        EXPECT_NEAR(c[n], 0.0, TOL) << "c[" << n << "] should be 0 for constant u_h";
}

// ─────────────────────────────────────────────────────────────────────────────
//  modalCoeffs: linear polynomial u_h(x) = x = P_1(x)
// ─────────────────────────────────────────────────────────────────────────────
TEST(ModalCoeffsTest, LinearPolynomial)
{
    DG::ReferenceElement<double> ref(4);
    DG::PadeLegendreSolver<double> solver(ref, 2, 2);

    const int Np = ref.numDOF();
    std::vector<double> nodal(Np);
    for (int i = 0; i < Np; ++i)
        nodal[i] = ref.nodes()[i];   // u(x) = x = P_1(x)

    const auto c = solver.modalCoeffs(nodal.data());

    EXPECT_NEAR(c[0], 0.0, TOL);
    EXPECT_NEAR(c[1], 1.0, TOL);
    for (int n = 2; n < Np; ++n)
        EXPECT_NEAR(c[n], 0.0, TOL) << "c[" << n << "] should be 0 for P_1";
}

// ─────────────────────────────────────────────────────────────────────────────
//  modalCoeffs: quadratic P_2(x) = (3x^2 - 1)/2
// ─────────────────────────────────────────────────────────────────────────────
TEST(ModalCoeffsTest, QuadraticLegendrePolynomial)
{
    DG::ReferenceElement<double> ref(4);
    DG::PadeLegendreSolver<double> solver(ref, 2, 2);

    const int Np = ref.numDOF();
    std::vector<double> nodal(Np);
    for (int i = 0; i < Np; ++i) {
        const double x = ref.nodes()[i];
        nodal[i] = (3.0 * x * x - 1.0) / 2.0;   // P_2(x)
    }

    const auto c = solver.modalCoeffs(nodal.data());

    EXPECT_NEAR(c[0], 0.0, TOL);
    EXPECT_NEAR(c[1], 0.0, TOL);
    EXPECT_NEAR(c[2], 1.0, TOL);
    for (int n = 3; n < Np; ++n)
        EXPECT_NEAR(c[n], 0.0, TOL) << "c[" << n << "] should be 0 for P_2";
}

// ─────────────────────────────────────────────────────────────────────────────
//  M = 0: Padé approximant reduces to the degree-N DG polynomial
//
//  With L = N, M = 0 the denominator is identically 1 and the numerator
//  equals the DG polynomial itself.  Evaluating at the GLL nodes must
//  reproduce the original nodal values.
// ─────────────────────────────────────────────────────────────────────────────
TEST(PadeLegendreSolverTest, PolyApproximationWithM0)
{
    const int N = 4;
    DG::ReferenceElement<double> ref(N);
    DG::PadeLegendreSolver<double> solver(ref, N, 0);

    const int Np = ref.numDOF();
    std::vector<double> nodal(Np);
    for (int i = 0; i < Np; ++i)
        nodal[i] = std::sin(ref.nodes()[i]);

    // Reconstruct at the GLL nodes; must reproduce nodal values
    std::vector<double> gll(Np);
    for (int i = 0; i < Np; ++i) gll[i] = ref.nodes()[i];
    const auto vals = solver.reconstruct(nodal.data(), gll);

    for (int i = 0; i < Np; ++i)
        EXPECT_NEAR(vals[i], nodal[i], TOL_LOOSE)
            << "M=0 reconstruction mismatch at node " << i;
}

// ─────────────────────────────────────────────────────────────────────────────
//  [1/1] Padé of 1/(1+0.5x) — analytical verification
//
//  The exact rational function is  1/(1 + 0.5x)  which has the [1/1]
//  Padé–Legendre representation  P_0 / (P_0 + 0.5 P_1)  =  1 / (1 + 0.5x).
//
//  Analytically (verified by hand):
//    q_1 = 0.5 ,  p_0 = 1  (exact regardless of polynomial approximation error).
// ─────────────────────────────────────────────────────────────────────────────
TEST(PadeLegendreSolverTest, Rational11PadeExactCoefficients)
{
    const int N = 2;
    DG::ReferenceElement<double> ref(N);
    DG::PadeLegendreSolver<double> solver(ref, 1, 1);

    const int Np = ref.numDOF();
    std::vector<double> nodal(Np);
    for (int i = 0; i < Np; ++i) {
        const double x = ref.nodes()[i];
        nodal[i] = 1.0 / (1.0 + 0.5 * x);
    }

    const auto c      = solver.modalCoeffs(nodal.data());
    const auto approx = solver.buildApproximant(c);

    EXPECT_NEAR(approx.q_coeffs[1], 0.5, TOL_PADE)
        << "Denominator coefficient q_1 should be 0.5";
    EXPECT_NEAR(approx.p_coeffs[0], 1.0, TOL_PADE)
        << "Numerator coefficient p_0 should be 1.0";
    EXPECT_TRUE(approx.valid) << "Approximant should be marked valid";
}

// ─────────────────────────────────────────────────────────────────────────────
//  [1/1] Padé of 1/(1+0.5x) — evaluate recovers the rational function
// ─────────────────────────────────────────────────────────────────────────────
TEST(PadeLegendreSolverTest, Rational11PadeEvaluationAccuracy)
{
    const int N = 2;
    DG::ReferenceElement<double> ref(N);
    DG::PadeLegendreSolver<double> solver(ref, 1, 1);

    const int Np = ref.numDOF();
    std::vector<double> nodal(Np);
    for (int i = 0; i < Np; ++i) {
        const double x = ref.nodes()[i];
        nodal[i] = 1.0 / (1.0 + 0.5 * x);
    }

    const auto c      = solver.modalCoeffs(nodal.data());
    const auto approx = solver.buildApproximant(c);

    // The Padé should reproduce 1/(1+0.5x) well at interior and boundary points
    for (double x : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
        const double exact = 1.0 / (1.0 + 0.5 * x);
        const double pade  = solver.evaluate(approx, x);
        EXPECT_NEAR(pade, exact, 1e-3)
            << "[1/1] Padé evaluation error too large at x=" << x;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Padé reduces L∞ error vs polynomial for 1/(1+0.5x)
//
//  With N = 4 and L = M = 2, the diagonal [2/2] Padé–Legendre approximant
//  of 1/(1+0.5x) should be more accurate than the degree-4 polynomial
//  approximation at interior test points.
// ─────────────────────────────────────────────────────────────────────────────
TEST(PadeLegendreSolverTest, PadeReducesErrorForRationalFunction)
{
    const int N = 4;
    DG::ReferenceElement<double> ref(N);
    DG::PadeLegendreSolver<double> solver(ref, 2, 2);

    const int Np = ref.numDOF();
    std::vector<double> nodal(Np);
    for (int i = 0; i < Np; ++i) {
        const double x = ref.nodes()[i];
        nodal[i] = 1.0 / (1.0 + 0.5 * x);
    }

    const auto c      = solver.modalCoeffs(nodal.data());
    const auto approx = solver.buildApproximant(c);

    // Sample at many interior test points (avoid the GLL nodes themselves)
    const int nTest = 50;
    double poly_err = 0.0, pade_err = 0.0;
    for (int i = 0; i < nTest; ++i) {
        const double x     = -0.99 + 1.98 * i / (nTest - 1);
        const double exact = 1.0 / (1.0 + 0.5 * x);

        // Polynomial evaluation: sum c_n P_n(x)
        double poly = 0.0;
        for (int n = 0; n < Np; ++n)
            poly += c[n] * DG::ReferenceElement<double>::legendreP(n, x);

        const double pade = solver.evaluate(approx, x);

        poly_err = std::max(poly_err, std::abs(poly - exact));
        pade_err = std::max(pade_err, std::abs(pade - exact));
    }

    EXPECT_LT(pade_err, poly_err)
        << "Padé [2/2] should have smaller L∞ error than degree-4 polynomial "
           "for 1/(1+0.5x)  (poly_err=" << poly_err
        << "  pade_err=" << pade_err << ")";
}

// ─────────────────────────────────────────────────────────────────────────────
//  reconstruct(FieldVector) returns same shape and finite values
// ─────────────────────────────────────────────────────────────────────────────
TEST(PadeLegendreSolverTest, ReconstructFieldVectorShape)
{
    const int K = 5, N = 4;
    DG::ReferenceElement<double> ref(N);
    DG::PadeLegendreSolver<double> solver(ref, 2, 2);

    DG::FieldVector<double> u(K, ref.numDOF());
    for (int k = 0; k < K; ++k)
        for (int i = 0; i < ref.numDOF(); ++i)
            u.elementPtr(k)[i] = std::sin(ref.nodes()[i]);

    const auto u_pade = solver.reconstruct(u);

    EXPECT_EQ(u_pade.numElements(), K);
    EXPECT_EQ(u_pade.numDOF(), ref.numDOF());

    for (int k = 0; k < K; ++k)
        for (int i = 0; i < ref.numDOF(); ++i)
            EXPECT_TRUE(std::isfinite(u_pade.elementPtr(k)[i]))
                << "Non-finite value at k=" << k << " i=" << i;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Constructor throws for invalid L, M
// ─────────────────────────────────────────────────────────────────────────────
TEST(PadeLegendreSolverTest, ConstructorThrowsForLplusMgtN)
{
    DG::ReferenceElement<double> ref(3);
    EXPECT_THROW(DG::PadeLegendreSolver<double>(ref, 2, 2),
                 std::invalid_argument)
        << "L+M=4 > N=3 should throw";
}

TEST(PadeLegendreSolverTest, ConstructorThrowsForNegativeM)
{
    DG::ReferenceElement<double> ref(4);
    EXPECT_THROW(DG::PadeLegendreSolver<double>(ref, 2, -1),
                 std::invalid_argument)
        << "Negative M should throw";
}

// ─────────────────────────────────────────────────────────────────────────────
//  Denominator is positive (no spurious poles) for a smooth function
// ─────────────────────────────────────────────────────────────────────────────
TEST(PadeLegendreSolverTest, DenominatorPositiveForSmoothFunction)
{
    const int N = 4;
    DG::ReferenceElement<double> ref(N);
    DG::PadeLegendreSolver<double> solver(ref, 2, 2);

    const int Np = ref.numDOF();
    std::vector<double> nodal(Np);
    for (int i = 0; i < Np; ++i)
        nodal[i] = std::exp(ref.nodes()[i]);  // smooth, no poles in [-1,1]

    const auto c      = solver.modalCoeffs(nodal.data());
    const auto approx = solver.buildApproximant(c);

    // Evaluate denominator Q at many points; must stay positive
    for (int i = 0; i < 50; ++i) {
        const double x = -1.0 + 2.0 * i / 49.0;
        double Q = approx.q_coeffs[0];
        for (int m = 1; m <= approx.M; ++m)
            Q += approx.q_coeffs[m]
                 * DG::ReferenceElement<double>::legendreP(m, x);
        EXPECT_GT(Q, 0.0) << "Denominator should be positive at x=" << x;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  writePadeVTK: output has the correct number of points (K * refineFactor * Np)
// ─────────────────────────────────────────────────────────────────────────────
TEST(WritePadeVTKTest, CorrectPointCount)
{
    const int K = 4, N = 4;
    const int refineFactor = 4;

    DG::Mesh<double>             mesh   = DG::Mesh<double>::uniform(0.0, 1.0, K);
    DG::ReferenceElement<double> ref(N);
    DG::PadeLegendreSolver<double> solver(ref, N / 2, N - N / 2);

    const int Np = ref.numDOF();
    DG::FieldVector<double> u(K, Np);
    for (int k = 0; k < K; ++k)
        for (int i = 0; i < Np; ++i)
            u.elementPtr(k)[i] = std::sin(ref.nodes()[i]);

    const std::string tmp = "/tmp/test_pade_vtk_points.vtk";
    ASSERT_NO_THROW(DG::writePadeVTK(mesh, ref, u, solver, tmp));

    // Count how many coordinate lines appear after "POINTS N double"
    std::ifstream fin(tmp);
    ASSERT_TRUE(fin.is_open()) << "VTK file was not created";

    const int expectedPts = K * refineFactor * Np;
    std::string line;
    int pointsFound = 0;
    bool inPoints = false;
    while (std::getline(fin, line)) {
        if (line.rfind("POINTS", 0) == 0) {
            // Extract declared count
            int declared = 0;
            std::istringstream ss(line);
            std::string tok;
            ss >> tok >> declared;
            EXPECT_EQ(declared, expectedPts) << "POINTS header count mismatch";
            inPoints = true;
            continue;
        }
        if (inPoints) {
            if (line.empty() || line.rfind("CELLS", 0) == 0) break;
            ++pointsFound;
        }
    }
    EXPECT_EQ(pointsFound, expectedPts) << "Actual point rows mismatch";
}

// ─────────────────────────────────────────────────────────────────────────────
//  writePadeVTK: output contains only finite values
// ─────────────────────────────────────────────────────────────────────────────
TEST(WritePadeVTKTest, AllValuesFinite)
{
    const int K = 3, N = 4;

    DG::Mesh<double>             mesh   = DG::Mesh<double>::uniform(-1.0, 1.0, K);
    DG::ReferenceElement<double> ref(N);
    DG::PadeLegendreSolver<double> solver(ref, 2, 2);

    const int Np = ref.numDOF();
    DG::FieldVector<double> u(K, Np);
    for (int k = 0; k < K; ++k)
        for (int i = 0; i < Np; ++i)
            u.elementPtr(k)[i] = std::exp(ref.nodes()[i]);

    const std::string tmp = "/tmp/test_pade_vtk_finite.vtk";
    ASSERT_NO_THROW(DG::writePadeVTK(mesh, ref, u, solver, tmp));

    // Scan POINT_DATA section and check every scalar value is finite
    std::ifstream fin(tmp);
    ASSERT_TRUE(fin.is_open());

    std::string line;
    bool inData = false;
    int scalarCount = 0;
    while (std::getline(fin, line)) {
        if (line.rfind("LOOKUP_TABLE", 0) == 0) { inData = true; continue; }
        if (!inData || line.empty()) continue;
        double val = 0.0;
        std::istringstream ss(line);
        ss >> val;
        EXPECT_TRUE(std::isfinite(val))
            << "Non-finite scalar at line: " << line;
        ++scalarCount;
    }
    EXPECT_GT(scalarCount, 0) << "No scalar data found in VTK file";
}

// ─────────────────────────────────────────────────────────────────────────────
//  writePadeVTK: throws when refineFactor < 1
// ─────────────────────────────────────────────────────────────────────────────
TEST(WritePadeVTKTest, ThrowsForZeroRefineFactor)
{
    DG::Mesh<double>             mesh   = DG::Mesh<double>::uniform(0.0, 1.0, 2);
    DG::ReferenceElement<double> ref(4);
    DG::PadeLegendreSolver<double> solver(ref, 2, 2);

    const int Np = ref.numDOF();
    DG::FieldVector<double> u(2, Np);

    EXPECT_THROW(DG::writePadeVTK(mesh, ref, u, solver,
                                   "/tmp/test_pade_bad.vtk", "u", 0.0, 0),
                 std::invalid_argument);
}
