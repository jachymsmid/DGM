// NumericalFlux.hpp
#pragma once

#include <cmath>
#include <functional>
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
  explicit UpwindFlux(std::function< Real( Real u ) > advection_speed, std::function< Real( Real u ) > physical_flux) : advection_speed_(std::move(advection_speed)), physical_flux_(std::move(physical_flux)) {}

  Real compute(Real u_minus, Real u_plus, Real n_outward) const override
  {
    Real C = std::max(advection_speed_(u_minus), advection_speed_(u_plus));
    return ( C * n_outward < 0) ? physical_flux_(u_minus) : physical_flux_(u_plus);
  }
  std::function< Real( Real u ) > advection_speed_;
  std::function< Real( Real u ) > physical_flux_;
};

// Lax-Friedrichs flux: f* = 0.5*(f(u-) + f(u+)) - 0.5*C*(u+ - u-)
template< class Real >
struct LaxFriedrichsFlux : NumericalFlux< Real >
{
  explicit LaxFriedrichsFlux( std::function< Real( Real u ) > advection_speed, std::function< Real( Real u ) > physical_flux ) : advection_speed_(std::move(advection_speed)), physical_flux_(std::move(physical_flux)) {}

  Real compute(Real u_minus, Real u_plus, Real n_outward) const override
  {
    Real C = std::max(std::fabs(advection_speed_(u_minus)), std::fabs(advection_speed_(u_plus)));
    return Real(0.5) * (physical_flux_(u_minus) + physical_flux_(u_plus)) + Real(0.5) * C * n_outward * (u_plus - u_minus);
  }
  std::function< Real( Real u ) > advection_speed_;
  std::function< Real( Real u ) > physical_flux_;
};

} // namespace DG
