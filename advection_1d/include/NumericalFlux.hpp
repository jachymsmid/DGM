/**
 * @file NumericalFlux.hpp
 * @brief Numerical flux interface and common flux implementations.
 *
 * Defines an abstract NumericalFlux interface and multiple concrete
 * implementations (UpwindFlux, LaxFriedrichsFlux, GodunovFlux, RoeFlux)
 * used to compute interface fluxes from interior/exterior states.
 */

#pragma once

#include <TNL/Math.h>
#include <functional>
namespace DG {

// ---------------------- Base numerical flux class ---------------------------
// Abstract base (crtp-free for simplicity; use std::function or templates)
template< class Real >
/**
 * @struct NumericalFlux
 * @brief Abstract interface for numerical flux computations at interfaces.
 *
 * Subclasses implement compute(u_minus,u_plus,n_outward) to return the
 * consistent inter-element flux for the DG surface term.
 */
struct NumericalFlux
{
  virtual ~NumericalFlux() = default;
  // f*: numerical flux at interface given the interior state u_minus (value
  // from the element on the "left" side of the interface), the exterior state
  // u_plus (right side), and the outward normal from the left element
  // (n_outward == -1 for the left face, +1 for the right face). The function
  // returns the scalar numerical flux f* (already oriented along the 1D
  // coordinate) to be used directly in the surface term assembly.
  virtual Real compute(Real u_minus, Real u_plus, Real n_outward) const = 0;
};

// -------------------------------- Upwind ------------------------------------
// f*(u) = f(u-)
template< class Real >
/**
 * @struct UpwindFlux
 * @brief Upwind numerical flux using the sign of the advection speed.
 *
 * Chooses the upwind state and evaluates the physical flux there.
 */
struct UpwindFlux : NumericalFlux< Real >
{
  // constructor
  explicit UpwindFlux(std::function< Real( Real u ) > advection_speed, std::function< Real( Real u ) > physical_flux) : advection_speed_(std::move(advection_speed)), physical_flux_(std::move(physical_flux)) {}

  Real compute(Real u_int, Real u_ext, Real n_outward) const override
  {
    // u_minus is the interior (left) state, u_plus is the exterior (right)
    // If a * n_outward >= 0 the characteristic leaves the interior cell so
    // the interior state determines the flux; otherwise use the exterior.
    Real C = advection_speed_(u_int);
    if (C * n_outward >= Real(0))
      return physical_flux_(u_int);
    else
      return physical_flux_(u_ext);
  }

  // data members
  std::function< Real( Real u ) > advection_speed_;
  std::function< Real( Real u ) > physical_flux_;
};


// ----------------------------- Lax-Friedrichs -------------------------------
// f* = 0.5*(f(u-) + f(u+)) - 0.5*C*(u+ - u-)
template< class Real >
/**
 * @struct LaxFriedrichsFlux
 * @brief Local Lax–Friedrichs (Rusanov) flux.
 *
 * Blends arithmetic flux with a jump penalization scaled by an estimate
 * of the maximum wave speed.
 */
struct LaxFriedrichsFlux : NumericalFlux< Real >
{
  explicit LaxFriedrichsFlux( std::function< Real( Real u ) > advection_speed, std::function< Real( Real u ) > physical_flux ) : advection_speed_(std::move(advection_speed)), physical_flux_(std::move(physical_flux)) {}

  Real compute(Real u_minus, Real u_plus, Real n_outward) const override
  {
    // alpha is an estimate of the maximum wave speed across the interface
    Real alpha = TNL::argAbsMax(advection_speed_(u_minus), advection_speed_(u_plus));
    return Real(0.5) * (physical_flux_(u_minus) + physical_flux_(u_plus))
           - Real(0.5) * alpha * n_outward * (u_plus - u_minus);
  }

  std::function< Real( Real u ) > advection_speed_;
  std::function< Real( Real u ) > physical_flux_;
};

// -------------------------------- Godunov ------------------------------------
template< class Real >
/**
 * @struct GodunovFlux
 * @brief Godunov flux (exact Riemann solver for scalar problems).
 */
struct GodunovFlux : NumericalFlux< Real >
{
  // constructor
  explicit GodunovFlux(std::function< Real( Real u ) > advection_speed, std::function< Real( Real u ) > physical_flux) : advection_speed_(std::move(advection_speed)), physical_flux_(std::move(physical_flux)) {}

  Real compute(Real u_minus, Real u_plus, Real n_outward) const override
  {
    return ( u_minus < u_plus ) ? n_outward * TNL::min(physical_flux_(u_minus), physical_flux_(u_plus)) : n_outward * TNL::max(physical_flux_(u_minus), physical_flux_(u_plus));
  }

  // data members
  std::function< Real( Real u ) > advection_speed_;
  std::function< Real( Real u ) > physical_flux_;
};

// -------------------------------- Roe ---------------------------------------
template< class Real >
/**
 * @struct RoeFlux
 * @brief Roe approximate Riemann solver (linearized flux).
 */
struct RoeFlux : NumericalFlux< Real >
{
  // constructor
  explicit RoeFlux(std::function< Real( Real u ) > advection_speed, std::function< Real( Real u ) > physical_flux) : advection_speed_(std::move(advection_speed)), physical_flux_(std::move(physical_flux)) {}

  Real compute(Real u_minus, Real u_plus, Real n_outward) const override
  {
    Real C = (advection_speed_(u_minus) + advection_speed_(u_plus))/2;
    return Real(0.5) * (physical_flux_(u_minus) + physical_flux_(u_plus)) - Real(0.5) * C * n_outward * (u_plus - u_minus);
  }

  // data members
  std::function< Real( Real u ) > advection_speed_;
  std::function< Real( Real u ) > physical_flux_;
};
} // namespace DG
