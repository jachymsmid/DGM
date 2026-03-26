// NumericalFlux.hpp
#pragma once

namespace DG {

// Abstract base (crtp-free for simplicity; use std::function or templates)
template< class Real >
struct NumericalFlux
{
  virtual ~NumericalFlux() = default;
  // f*: numerical flux at interface given left state u_minus, right state u_plus,
  // and outward normal from the left element (+- 1 in 1D).
  virtual Real compute(Real u_minus, Real u_plus, Real n_outward) const = 0;
};

// Upwind flux for scalar advection: f*(u) = f(u-)
template< class Real >
struct UpwindFlux : NumericalFlux< Real >
{
  // constructor
  explicit UpwindFlux(Real a) : a_(a) {}

  Real compute(Real u_minus, Real u_plus, Real n_outward) const override
  {
    Real an = a_ * n_outward;
    return (an < 0) ? a_ * u_minus : a_ * u_plus;
  }
  Real a_;
};

// Lax-Friedrichs flux: f* = 0.5*(f(u-) + f(u+)) - 0.5*C*(u+ - u-)
template< class Real >
struct LaxFriedrichsFlux : NumericalFlux< Real >
{
  LaxFriedrichsFlux(Real a, Real C) : a_(a), C_(C) {}
  Real compute(Real u_minus, Real u_plus, Real n_outward) const override
  {
    return Real(0.5) * (a_*u_minus + a_*u_plus) - Real(0.5) * C_ * (u_plus - u_minus);
  }
  Real a_, C_;
};

} // namespace DG
