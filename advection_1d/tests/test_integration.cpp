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
