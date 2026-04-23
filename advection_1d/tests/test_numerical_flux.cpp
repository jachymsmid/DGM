/**
 * @file test_numerical_flux.cpp
 * @brief Unit tests for numerical flux implementations.
 *
 * Covers UpwindFlux, LaxFriedrichsFlux, GodunovFlux, and RoeFlux with
 * consistency checks and specific known-value cases.
 */

#include <gtest/gtest.h>
#include "NumericalFlux.hpp"
#include <functional>
#include <cmath>

static constexpr double TOL = 1e-12;

// ═══════════════════════════════════════════════════════════════════════════════
//  UpwindFlux
// ═══════════════════════════════════════════════════════════════════════════════

class UpwindFluxTest : public ::testing::Test
{
protected:
    // a = +1: characteristics travel right; f(u) = u
    std::function<double(double)> speedPos  = [](double /*u*/) { return  1.0; };
    std::function<double(double)> fluxPos_f = [](double u)     { return  u;   };
    DG::UpwindFlux<double> fluxPos{ speedPos, fluxPos_f };

    // a = -1: characteristics travel left; f(u) = -u
    std::function<double(double)> speedNeg  = [](double /*u*/) { return -1.0; };
    std::function<double(double)> fluxNeg_f = [](double u)     { return -u;   };
    DG::UpwindFlux<double> fluxNeg{ speedNeg, fluxNeg_f };
};

// Positive speed, right face (n=+1): C*n = 1*1 > 0 → use interior state
TEST_F(UpwindFluxTest, PositiveSpeedRightFaceUsesInterior)
{
    double u_int = 2.0, u_ext = 3.0, n = +1.0;
    // f(u_int) = 1 * 2 = 2
    EXPECT_NEAR(fluxPos.compute(u_int, u_ext, n), 2.0, TOL);
}

// Positive speed, left face (n=-1): C*n = 1*(-1) < 0 → use exterior state
TEST_F(UpwindFluxTest, PositiveSpeedLeftFaceUsesExterior)
{
    double u_int = 2.0, u_ext = 3.0, n = -1.0;
    // f(u_ext) = 1 * 3 = 3
    EXPECT_NEAR(fluxPos.compute(u_int, u_ext, n), 3.0, TOL);
}

// Negative speed, right face (n=+1): C*n = (-1)*1 < 0 → use exterior state
TEST_F(UpwindFluxTest, NegativeSpeedRightFaceUsesExterior)
{
    double u_int = 2.0, u_ext = 3.0, n = +1.0;
    // f(u_ext) = -1 * 3 = -3
    EXPECT_NEAR(fluxNeg.compute(u_int, u_ext, n), -3.0, TOL);
}

// Negative speed, left face (n=-1): C*n = (-1)*(-1) = 1 ≥ 0 → use interior
TEST_F(UpwindFluxTest, NegativeSpeedLeftFaceUsesInterior)
{
    double u_int = 2.0, u_ext = 3.0, n = -1.0;
    // f(u_int) = -1 * 2 = -2
    EXPECT_NEAR(fluxNeg.compute(u_int, u_ext, n), -2.0, TOL);
}

// Consistency: equal states give f(u) regardless of normal direction
TEST_F(UpwindFluxTest, ConsistencyRightFace)
{
    double u = 4.0;
    EXPECT_NEAR(fluxPos.compute(u, u, +1.0), 4.0, TOL);
    EXPECT_NEAR(fluxPos.compute(u, u, -1.0), 4.0, TOL);
}

TEST_F(UpwindFluxTest, ConsistencyNegativeSpeed)
{
    double u = -2.5;
    // f(u) = -(-2.5) = 2.5
    EXPECT_NEAR(fluxNeg.compute(u, u, +1.0), 2.5, TOL);
    EXPECT_NEAR(fluxNeg.compute(u, u, -1.0), 2.5, TOL);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  LaxFriedrichsFlux
// ═══════════════════════════════════════════════════════════════════════════════

class LaxFriedrichsFluxTest : public ::testing::Test
{
protected:
    // Burgers: f(u) = u^2/2, speed(u) = u
    std::function<double(double)> burgSpeed = [](double u) { return u; };
    std::function<double(double)> burgFlux  = [](double u) { return 0.5 * u * u; };
    DG::LaxFriedrichsFlux<double> flux{ burgSpeed, burgFlux };
};

// Consistency: u_minus = u_plus → dissipation vanishes → flux = f(u)
TEST_F(LaxFriedrichsFluxTest, ConsistencyEqualStates)
{
    double u = 3.0;
    // f(u) = u^2/2 = 4.5
    EXPECT_NEAR(flux.compute(u, u, +1.0), 4.5, TOL);
    EXPECT_NEAR(flux.compute(u, u, -1.0), 4.5, TOL);
}

TEST_F(LaxFriedrichsFluxTest, ConsistencyZeroState)
{
    EXPECT_NEAR(flux.compute(0.0, 0.0, +1.0), 0.0, TOL);
}

// For Burgers, u_minus=0, u_plus=2, n=+1:
//   C = argAbsMax(speed(0), speed(2)) = argAbsMax(0, 2) = 2
//   flux = 0.5*(f(0)+f(2)) + 0.5*C*n*(u_minus - u_plus)
//        = 0.5*(0+2) + 0.5*2*1*(0-2) = 1 - 2 = -1
TEST_F(LaxFriedrichsFluxTest, KnownValueBurgers)
{
    double C = 2.0; // |speed(u_plus)| = |2| = 2 is the larger
    double expected = 0.5 * (0.0 + 2.0) + 0.5 * C * 1.0 * (0.0 - 2.0);
    EXPECT_NEAR(flux.compute(0.0, 2.0, +1.0), expected, TOL);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  GodunovFlux
// ═══════════════════════════════════════════════════════════════════════════════

class GodunovFluxTest : public ::testing::Test
{
protected:
    // Linear advection: f(u)=u, speed=1
    std::function<double(double)> linSpeed  = [](double /*u*/) { return 1.0; };
    std::function<double(double)> linFlux_f = [](double u)     { return u;   };
    DG::GodunovFlux<double> fluxLinear{ linSpeed, linFlux_f };

    // Burgers: f(u)=u^2/2, speed=u
    std::function<double(double)> burgSpeed = [](double u) { return u; };
    std::function<double(double)> burgFlux  = [](double u) { return 0.5 * u * u; };
    DG::GodunovFlux<double> fluxBurgers{ burgSpeed, burgFlux };
};

// Consistency: equal states give f(u)
TEST_F(GodunovFluxTest, ConsistencyLinearAdvection)
{
    double u = 2.5;
    EXPECT_NEAR(fluxLinear.compute(u, u, +1.0), 2.5, TOL);
    EXPECT_NEAR(fluxLinear.compute(u, u, -1.0), 2.5, TOL);
}

// Linear advection, n=+1, u_minus < u_plus: n*u_minus < n*u_plus → min(f(u-),f(u+))
// f(u)=u, u_minus=1 < u_plus=3: min(1,3) = 1
TEST_F(GodunovFluxTest, LinearRarefactionRightFace)
{
    EXPECT_NEAR(fluxLinear.compute(1.0, 3.0, +1.0), 1.0, TOL);
}

// Linear advection, n=+1, u_minus > u_plus: max(f(u-),f(u+))
// u_minus=3 > u_plus=1: max(3,1) = 3
TEST_F(GodunovFluxTest, LinearShockRightFace)
{
    EXPECT_NEAR(fluxLinear.compute(3.0, 1.0, +1.0), 3.0, TOL);
}

// Burgers rarefaction: n=+1, u_minus=1 < u_plus=2
//   n*u_minus=1 < n*u_plus=2 → min(f(1),f(2)) = min(0.5, 2.0) = 0.5
TEST_F(GodunovFluxTest, BurgersRarefactionRightFace)
{
    EXPECT_NEAR(fluxBurgers.compute(1.0, 2.0, +1.0), 0.5, TOL);
}

// Burgers shock: n=+1, u_minus=2 > u_plus=1
//   n*u_minus=2 > n*u_plus=1 → max(f(2),f(1)) = max(2.0, 0.5) = 2.0
TEST_F(GodunovFluxTest, BurgersShockRightFace)
{
    EXPECT_NEAR(fluxBurgers.compute(2.0, 1.0, +1.0), 2.0, TOL);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  RoeFlux
// ═══════════════════════════════════════════════════════════════════════════════

class RoeFluxTest : public ::testing::Test
{
protected:
    // Linear advection: f(u)=u, speed=1
    std::function<double(double)> linSpeed  = [](double /*u*/) { return 1.0; };
    std::function<double(double)> linFlux_f = [](double u)     { return u;   };
    DG::RoeFlux<double> fluxLinear{ linSpeed, linFlux_f };

    // Burgers: f(u)=u^2/2, speed=u
    std::function<double(double)> burgSpeed = [](double u) { return u; };
    std::function<double(double)> burgFlux  = [](double u) { return 0.5 * u * u; };
    DG::RoeFlux<double> fluxBurgers{ burgSpeed, burgFlux };
};

// Consistency: equal states give f(u)
TEST_F(RoeFluxTest, ConsistencyLinearAdvection)
{
    double u = 3.0;
    EXPECT_NEAR(fluxLinear.compute(u, u, +1.0), 3.0, TOL);
    EXPECT_NEAR(fluxLinear.compute(u, u, -1.0), 3.0, TOL);
}

TEST_F(RoeFluxTest, ConsistencyBurgers)
{
    double u = 2.0;
    // f(u) = u^2/2 = 2.0
    EXPECT_NEAR(fluxBurgers.compute(u, u, +1.0), 2.0, TOL);
}

// Burgers with Roe flux: u_minus=1, u_plus=3, n=+1
//   alpha = (|speed(1)|+|speed(3)|)/2 = (1+3)/2 = 2
//   flux = 0.5*(f(1)+f(3)) + 0.5*alpha*n*(u_minus - u_plus)
//        = 0.5*(0.5+4.5) + 0.5*2*1*(1-3) = 2.5 - 2 = 0.5
TEST_F(RoeFluxTest, BurgersKnownValueRightFace)
{
    double u_minus = 1.0, u_plus = 3.0, n = +1.0;
    double alpha = (std::abs(u_minus) + std::abs(u_plus)) / 2.0;
    double expected = 0.5 * (0.5 + 4.5) + 0.5 * alpha * n * (u_minus - u_plus);
    EXPECT_NEAR(fluxBurgers.compute(u_minus, u_plus, n), expected, TOL);
}

// Consistency: flux is symmetric when states are equal
TEST_F(RoeFluxTest, SymmetryEqualStates)
{
    double u = 1.5;
    EXPECT_NEAR(fluxLinear.compute(u, u, +1.0),
                fluxLinear.compute(u, u, -1.0), TOL);
}
