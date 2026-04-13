// NumericalFlux.hpp
#pragma once

#include <TNL/Math.h>
#include <functional>
namespace DG {

// ---------------------- Base numerical flux class ---------------------------
// Abstract base (crtp-free for simplicity; use std::function or templates)
template< class Real >
struct NumericalFlux
{
  virtual ~NumericalFlux() = default;
  // f*: numerical flux at interface given left state u_minus, right state u_plus,
  // and outward normal from the left element (+- 1 in 1D).
  virtual Real compute(Real u_minus, Real u_plus, Real n_outward) const = 0;
};

// -------------------------------- Upwind ------------------------------------
// f*(u) = f(u-)
template< class Real >
struct UpwindFlux : NumericalFlux< Real >
{
  // constructor
  explicit UpwindFlux(std::function< Real( Real u ) > advection_speed, std::function< Real( Real u ) > physical_flux) : advection_speed_(std::move(advection_speed)), physical_flux_(std::move(physical_flux)) {}

  Real compute(Real u_minus, Real u_plus, Real n_outward) const override
  {
    Real C = advection_speed_(u_minus);
    return ( C * n_outward < 0) ? physical_flux_(u_minus) : physical_flux_(u_plus);
  }

  // data members
  std::function< Real( Real u ) > advection_speed_;
  std::function< Real( Real u ) > physical_flux_;
};


// ----------------------------- Lax-Friedrichs -------------------------------
// f* = 0.5*(f(u-) + f(u+)) - 0.5*C*(u+ - u-)
template< class Real >
struct LaxFriedrichsFlux : NumericalFlux< Real >
{
  explicit LaxFriedrichsFlux( std::function< Real( Real u ) > advection_speed, std::function< Real( Real u ) > physical_flux ) : advection_speed_(std::move(advection_speed)), physical_flux_(std::move(physical_flux)) {}

  Real compute(Real u_minus, Real u_plus, Real n_outward) const override
  {
    Real C = TNL::argAbsMax(advection_speed_(u_minus), advection_speed_(u_plus));
    return n_outward * Real(0.5) * (physical_flux_(u_minus) + physical_flux_(u_plus)) - Real(0.5) * C * (u_plus - u_minus);
  }
  std::function< Real( Real u ) > advection_speed_;
  std::function< Real( Real u ) > physical_flux_;
};

// -------------------------------- Godunov ------------------------------------
template< class Real >
struct GodunovFlux : NumericalFlux< Real >
{
  // constructor
  explicit GodunovFlux(std::function< Real( Real u ) > advection_speed, std::function< Real( Real u ) > physical_flux) : advection_speed_(std::move(advection_speed)), physical_flux_(std::move(physical_flux)) {}

  Real compute(Real u_minus, Real u_plus, Real n_outward) const override
  {
    return ( u_minus < u_plus ) ? TNL::min(physical_flux_(u_minus), physical_flux_(u_plus)) : TNL::max(physical_flux_(u_minus), physical_flux_(u_plus));
  }

  // data members
  std::function< Real( Real u ) > advection_speed_;
  std::function< Real( Real u ) > physical_flux_;
};


} // namespace DG
