#include <gtest/gtest.h>
#include "FieldVector.hpp"
#include "Integrator.hpp"
#include <cmath>

// tolerance
double TOL = 1e-4;

// test on the problem y' = -y, y(0) = 1
TEST(ERKTest, ExponentialDecay)
{
  double t_end = 1.0;
  double dt = 0.01;
  double y_0 = 1.0;

  DG::FieldVector<double> y(1, 1);

  // y(0) = 1
  y.elementPtr(0)[0] = y_0;

  auto rhs = [&](const DG::FieldVector<double>& u_in,
                       DG::FieldVector<double>& u_out,
                 const double &time)
  {
    u_out.elementPtr(0)[0] = -u_in.elementPtr(0)[0];
  };

  DG::ERK<double> rk(rhs, 1, 1);
  rk.integrate(y, 0.0, t_end, dt);

  EXPECT_NEAR(y.elementPtr(0)[0], std::exp(-1), TOL);
}

// test the low storage explicit runge-kutta
TEST(LSERKTest, ExponentialDecay)
{
  double t_end = 1.0;
  double dt = 0.1;
  double y_0 = 1.0;

  DG::FieldVector<double> y(1, 1);

  // y(0) = 1
  y.elementPtr(0)[0] = y_0;

  auto rhs = [&](const DG::FieldVector<double>& u_in,
                       DG::FieldVector<double>& u_out,
                 const double &time)
  {
    u_out.elementPtr(0)[0] = -u_in.elementPtr(0)[0];
  };

  DG::LSERK<double> rk(rhs, 1, 1);
  rk.integrate(y, 0.0, t_end, dt);

  EXPECT_NEAR(y.elementPtr(0)[0], std::exp(-1), TOL);
}

// test on the problem y' = -y, y(0) = 1
TEST(SSPRKTest, ExponentialDecay)
{
  double t_end = 1.0;
  double dt = 0.01;
  double y_0 = 1.0;

  DG::FieldVector<double> y(1, 1);

  // y(0) = 1
  y.elementPtr(0)[0] = y_0;

  auto rhs = [&](const DG::FieldVector<double>& u_in,
                       DG::FieldVector<double>& u_out,
                 const double &time)
  {
    u_out.elementPtr(0)[0] = -u_in.elementPtr(0)[0];
  };

  DG::SSPRK<double> rk(rhs, 1, 1);
  rk.integrate(y, 0.0, t_end, dt);

  EXPECT_NEAR(y.elementPtr(0)[0], std::exp(-1), TOL);
}

// ── Fixture helpers ────────────────────────────────────────────────────────────
// Shared RHS for y' = -y  used by error-handling and ancillary tests.
static auto makeDecayRHS()
{
    return [](const DG::FieldVector<double>& u_in,
                    DG::FieldVector<double>& u_out,
              const double& /*time*/)
    {
        u_out.elementPtr(0)[0] = -u_in.elementPtr(0)[0];
    };
}

// ── Error handling: step with non-positive dt ──────────────────────────────────
TEST(ERKErrorTest, StepNonPositiveDtThrows)
{
    auto rhs = makeDecayRHS();
    DG::ERK<double> rk(rhs, 1, 1);
    DG::FieldVector<double> y(1, 1);
    y.elementPtr(0)[0] = 1.0;
    EXPECT_THROW(rk.step(y, 0.0,  0.0), std::invalid_argument);
    EXPECT_THROW(rk.step(y, -0.1, 0.0), std::invalid_argument);
}

TEST(LSERKErrorTest, StepNonPositiveDtThrows)
{
    auto rhs = makeDecayRHS();
    DG::LSERK<double> rk(rhs, 1, 1);
    DG::FieldVector<double> y(1, 1);
    y.elementPtr(0)[0] = 1.0;
    EXPECT_THROW(rk.step(y, 0.0,  0.0), std::invalid_argument);
    EXPECT_THROW(rk.step(y, -0.1, 0.0), std::invalid_argument);
}

TEST(SSPRKErrorTest, StepNonPositiveDtThrows)
{
    auto rhs = makeDecayRHS();
    DG::SSPRK<double> rk(rhs, 1, 1);
    DG::FieldVector<double> y(1, 1);
    y.elementPtr(0)[0] = 1.0;
    EXPECT_THROW(rk.step(y, 0.0,  0.0), std::invalid_argument);
    EXPECT_THROW(rk.step(y, -0.1, 0.0), std::invalid_argument);
}

// ── Error handling: integrate with invalid time range ────────────────────────
TEST(ERKErrorTest, IntegrateInvalidRangeThrows)
{
    auto rhs = makeDecayRHS();
    DG::ERK<double> rk(rhs, 1, 1);
    DG::FieldVector<double> y(1, 1);
    EXPECT_THROW(rk.integrate(y, 1.0, 0.5, 0.1), std::invalid_argument);
    EXPECT_THROW(rk.integrate(y, 1.0, 1.0, 0.1), std::invalid_argument);
}

TEST(LSERKErrorTest, IntegrateInvalidRangeThrows)
{
    auto rhs = makeDecayRHS();
    DG::LSERK<double> rk(rhs, 1, 1);
    DG::FieldVector<double> y(1, 1);
    EXPECT_THROW(rk.integrate(y, 1.0, 0.5, 0.1), std::invalid_argument);
    EXPECT_THROW(rk.integrate(y, 1.0, 1.0, 0.1), std::invalid_argument);
}

TEST(SSPRKErrorTest, IntegrateInvalidRangeThrows)
{
    auto rhs = makeDecayRHS();
    DG::SSPRK<double> rk(rhs, 1, 1);
    DG::FieldVector<double> y(1, 1);
    EXPECT_THROW(rk.integrate(y, 1.0, 0.5, 0.1), std::invalid_argument);
    EXPECT_THROW(rk.integrate(y, 1.0, 1.0, 0.1), std::invalid_argument);
}

// ── computeDt: zero wave speed throws ────────────────────────────────────────
TEST(ERKComputeDtTest, ZeroWaveSpeedThrows)
{
    EXPECT_THROW(DG::ERK<double>::computeDt(0.1, 0.0, 2), std::invalid_argument);
}

TEST(LSERKComputeDtTest, ZeroWaveSpeedThrows)
{
    EXPECT_THROW(DG::LSERK<double>::computeDt(0.1, 0.0, 2), std::invalid_argument);
}

TEST(SSPRKComputeDtTest, ZeroWaveSpeedThrows)
{
    EXPECT_THROW(DG::SSPRK<double>::computeDt(0.1, 0.0, 2), std::invalid_argument);
}

// ── computeDt: known value ─────────────────────────────────────────────────────
// dt = cfl * x_min / (|a| * (2N+1))
// With x_min=1.0, a=1.0, N=1, cfl=0.4: dt = 0.4 / 3
TEST(ERKComputeDtTest, KnownValue)
{
    double dt = DG::ERK<double>::computeDt(1.0, 1.0, 1, 0.4);
    EXPECT_NEAR(dt, 0.4 / 3.0, 1e-12);
}

TEST(SSPRKComputeDtTest, KnownValue)
{
    // dt = 0.4 * 0.5 / (2.0 * 5) = 0.02
    double dt = DG::SSPRK<double>::computeDt(0.5, 2.0, 2, 0.4);
    EXPECT_NEAR(dt, 0.02, 1e-12);
}

// ── numPoints ─────────────────────────────────────────────────────────────────
TEST(ERKTest, NumPoints)
{
    auto rhs = makeDecayRHS();
    DG::ERK<double> rk(rhs, 3, 4);
    EXPECT_EQ(rk.numPoints(), 12);
}

TEST(LSERKTest, NumPoints)
{
    auto rhs = makeDecayRHS();
    DG::LSERK<double> rk(rhs, 5, 2);
    EXPECT_EQ(rk.numPoints(), 10);
}

TEST(SSPRKTest, NumPoints)
{
    auto rhs = makeDecayRHS();
    DG::SSPRK<double> rk(rhs, 2, 6);
    EXPECT_EQ(rk.numPoints(), 12);
}

// ── Callback: invoked once per step ──────────────────────────────────────────
TEST(ERKTest, CallbackInvokedEachStep)
{
    auto rhs = makeDecayRHS();
    int callCount = 0;
    auto cb = [&](double /*t*/,
                  const DG::FieldVector<double>& /*u*/,
                  int /*step*/) -> bool
    {
        ++callCount;
        return true;
    };

    DG::FieldVector<double> y(1, 1);
    y.elementPtr(0)[0] = 1.0;
    DG::ERK<double> rk(rhs, 1, 1, cb);
    rk.integrate(y, 0.0, 1.0, 0.1);

    EXPECT_EQ(callCount, 10);
}

// ── Callback: returning false stops integration early ────────────────────────
TEST(ERKTest, CallbackCanStopEarly)
{
    auto rhs = makeDecayRHS();
    int callCount = 0;
    const int stopAfter = 3;

    auto cb = [&](double /*t*/,
                  const DG::FieldVector<double>& /*u*/,
                  int step) -> bool
    {
        ++callCount;
        return step < stopAfter;
    };

    DG::FieldVector<double> y(1, 1);
    y.elementPtr(0)[0] = 1.0;
    DG::ERK<double> rk(rhs, 1, 1, cb);
    rk.integrate(y, 0.0, 1.0, 0.1);  // would be 10 steps without early stop

    EXPECT_EQ(callCount, stopAfter);
}

TEST(LSERKTest, CallbackCanStopEarly)
{
    auto rhs = makeDecayRHS();
    int callCount = 0;
    const int stopAfter = 2;

    auto cb = [&](double /*t*/,
                  const DG::FieldVector<double>& /*u*/,
                  int step) -> bool
    {
        ++callCount;
        return step < stopAfter;
    };

    DG::FieldVector<double> y(1, 1);
    y.elementPtr(0)[0] = 1.0;
    DG::LSERK<double> rk(rhs, 1, 1, cb);
    rk.integrate(y, 0.0, 1.0, 0.1);

    EXPECT_EQ(callCount, stopAfter);
}

TEST(SSPRKTest, CallbackCanStopEarly)
{
    auto rhs = makeDecayRHS();
    int callCount = 0;
    const int stopAfter = 4;

    auto cb = [&](double /*t*/,
                  const DG::FieldVector<double>& /*u*/,
                  int step) -> bool
    {
        ++callCount;
        return step < stopAfter;
    };

    DG::FieldVector<double> y(1, 1);
    y.elementPtr(0)[0] = 1.0;
    DG::SSPRK<double> rk(rhs, 1, 1, cb);
    rk.integrate(y, 0.0, 1.0, 0.1);

    EXPECT_EQ(callCount, stopAfter);
}
