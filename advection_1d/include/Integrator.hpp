/**
 * @file Integrator.hpp
 * @brief Time integrators for advancing the semi-discrete DG system.
 *
 * Declares an abstract Integrator interface and several concrete explicit
 * Runge--Kutta schemes used by the solver (ERK, LSERK, SSPRK). Integrators
 * accept an RHS callback and operate on FieldVector instances.
 */

#pragma once

#include "FieldVector.hpp"
#include <TNL/Math.h>
#include <algorithm>
#include <array>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace DG {

//  callback type: called after every accepted step with (t, u, step_index).
template<typename Real, typename Device, typename Index>
using StepCallback = std::function< bool(Real t, const FieldVector<Real,Device,Index>& u, Index step) >;

// generic integrator class
template
<
  typename Real   = double,
  typename Device = TNL::Devices::Host,
  typename Index  = int
>
/**
 * @class Integrator
 * @brief Abstract base for time integrators advancing the semi-discrete DG system.
 *
 * Subclasses implement step() and integrate(). Integrators operate on
 * FieldVector instances and accept an RHS callback describing the spatial
 * operator.
 */
class Integrator
{
public:
    using Field = FieldVector<Real, Device, Index>;
    using RHSFunc = std::function<void(const Field&, Field&, const Real&)>;
    using Callback = StepCallback<Real, Device, Index>;

    virtual ~Integrator() = default;

    // single step
    /**
     * @brief Perform a single time step updating u in-place.
     * @param u solution vector to be advanced
     * @param dt time step size
     * @param t_in current time
     */
    virtual void step(Field& u, Real dt, Real t_in) = 0;

    // integrate, multiple steps
    /**
     * @brief Integrate from t0 to t_end using fixed step dt (last step clipped).
     * @return number of steps performed
     */
    virtual Index integrate(Field& u, Real t0, Real t_end, Real dt) = 0;

    // getters
    /**
     * @brief Number of steps performed by the integrator.
     */
    virtual Index numSteps() const = 0;
    /**
     * @brief Total number of DOFs managed by the integrator (K * Np).
     */
    virtual Index numPoints() const = 0;
    /**
     * @brief Current simulation time tracked by the integrator.
     */
    virtual Real currentTime() const = 0;

protected:
    // scale a field by a scalar
    void scaleField_(Field& f, Real s) const
    {
        const Index total = numPoints();
        Real* data = f.data().getData();
        for (Index i = 0; i < total; ++i)
        {
          data[i] *= s;
        }
    }

    // res = base + alpha * delta
    void addScaled_(const Field& base,
                    const Field& delta,
                    Real alpha,
                    Field& out) const
    {
      const Index total = numPoints();
      const Real* b = base.data().getData();
      const Real* d = delta.data().getData();
      Real* o = out.data().getData();
      for (Index i = 0; i < total; ++i)
      {
        o[i] = b[i] + alpha * d[i];
      }
    }
};


// ------------------------------- ERK4 ---------------------------------------
// simple explicit four stage Runge-Kutta
template<typename Real   = double,
         typename Device = TNL::Devices::Host,
         typename Index  = int>
/**
 * @class ERK
 * @brief Classic explicit four-stage Runge–Kutta integrator (ERK4).
 *
 * Implements a straightforward RK4 integration using temporary stage
 * storage sized by K * Np.
 */
class ERK : Integrator< Real, Device, Index > // inherit protected methods
{
public:
    using Field = FieldVector<Real, Device, Index>;
    using RHSFunc = std::function<void(const Field&, Field&, const Real&)>;
    using Callback = StepCallback<Real, Device, Index>;

    // op       – spatial residual (must outlive this integrator)
    // K, Np    – element count and nodes-per-element (needed to size scratch)
    // callback – optional: called after every step, return false to stop early
    /**
     * @brief Construct ERK integrator.
     * @param rhs Spatial residual callback: (u, out, time)
     * @param K number of elements
     * @param Np nodes per element
     * @param callback optional per-step callback; return false to stop
     */
    explicit ERK(RHSFunc rhs,
                 Index K, Index Np,
                 Callback callback = nullptr)
        : rhs_(std::move(rhs)),
          K_(K), Np_(Np),
          callback_(std::move(callback)),
          k1_(K, Np), k2_(K, Np), k3_(K, Np), k4_(K, Np),
          tmp_(K, Np) {}

    // single step
    void step(Field& u, Real dt, Real t_in) override
    {
      if (dt <= Real(0))
        throw std::invalid_argument("DG::Integrator::step: dt must be > 0");

      // k1 = L(u^n, t_n)
      rhs_(u, k1_, t_in + 0.5 * dt);

      // k2 = L(u^n + 1/2 dt k1, t_n + 1/2 dt)
      // tmp = u + 0.5*dt*k1
      this -> addScaled_(u, k1_, 0.5 * dt, tmp_);
      rhs_(tmp_, k2_, t_in + 0.5 * dt);

      // k3 = L(u^n + 1/2 dt k2, t_n + 1/2 dt)
      // tmp = u + 0.5*dt*k2
      this -> addScaled_(u, k2_, 0.5 * dt, tmp_);
      rhs_(tmp_, k3_, t_in + 0.5 * dt);

      // k4 = L(u^n + dt k3, t_n + dt)
      // tmp = u + dt * k3
      this -> addScaled_(u, k3_, dt, tmp_);
      rhs_(tmp_, k4_, t_in + dt );

      // combine: u += (k1 + 2*k2 + 2*k3 + k4) / 6
      combine_(u, k1_, k2_, k3_, k4_, dt);
    }

    // integrate: perform multiple steps
    Index integrate(Field& u, Real t0, Real t_end, Real dt) override
    {
      if (t_end <= t0)
        throw std::invalid_argument("Integrator::integrate: t_end must be > t0");

      Real t = t0;
      Index nstep = 0;

      while (t < t_end - Real(1e-12))
      {
        // clip the last step to reach t_end exactly
        Real dt_actual = std::min(dt, t_end - t);

        step(u, dt_actual, t);
        t += dt_actual;
        nstep += 1;

        // Invoke optional callback; stop early if it returns false
        if (callback_ && !callback_(t, u, nstep))
            break;
      }

      return nstep;
    }

    // max dt from the Hesthaven book
    static Real computeDt(Real x_min, Real max_wave_speed, int poly_order, Real cfl = Real(0.4))
    {
        if (TNL::abs(max_wave_speed) <= Real(0))
            throw std::invalid_argument("max_wave_speed must be positive");
        return cfl * x_min / (TNL::abs(max_wave_speed) * ( 2 * poly_order + 1 ));
    }

    // getters
    Index numSteps() const override { return lastStepCount_; }
    Real currentTime() const override { return currentTime_; }
    Index numPoints() const override { return K_ * Np_; }

private:
    // combine all the stages
    void combine_(Field& u,
                  const Field& k1,
                  const Field& k2,
                  const Field& k3,
                  const Field& k4,
                  const Real& dt) const
    {
      const Index total = K_ * Np_;
      Real* pu = u.data().getData();
      const Real* p1 = k1.data().getData();
      const Real* p2 = k2.data().getData();
      const Real* p3 = k3.data().getData();
      const Real* p4 = k4.data().getData();
      const Real sixth = Real(1) / Real(6);
      for (Index i = 0; i < total; ++i)
      {
        pu[i] += dt * sixth * (p1[i] + Real(2)*p2[i] + Real(2)*p3[i] + p4[i]);
      }
    }

    // RHS operator (function)
    RHSFunc rhs_;

    // number of cells and nodes
    Index K_, Np_;

    // the four stages
    Field k1_, k2_, k3_, k4_;
    Field tmp_;

    // callback function
    Callback callback_;

    // diagnostics variables, updated by step()
    Index lastStepCount_{ 0 };
    Real currentTime_{ 0 };
};

// ------------------------------- LSERK ---------------------------------------
// low storage explicit five stage Runge-Kutta
// see: https://ntrs.nasa.gov/api/citations/19940028444/downloads/19940028444.pdf
template<typename Real   = double,
         typename Device = TNL::Devices::Host,
         typename Index  = int>
class LSERK : Integrator< Real, Device, Index > // inherit protected methods
{
public:
    using Field = FieldVector<Real, Device, Index>;
    using RHSFunc = std::function<void(const Field&, Field&, const Real&)>;
    using Callback = StepCallback<Real, Device, Index>;

    // op       – spatial residual (must outlive this integrator)
    // K, Np    – element count and nodes-per-element (needed to size scratch)
    // callback – optional: called after every step, return false to stop early
    explicit LSERK(RHSFunc rhs,
                           Index K, Index Np,
                           Callback callback = nullptr)
        : rhs_(std::move(rhs))
        , K_(K), Np_(Np)
        , callback_(std::move(callback))
        , du_(K, Np)
        , tmp_(K, Np) {}

    // single step
    void step(Field& u, Real dt, Real t_in) override
    {
      if (dt <= Real(0))
        throw std::invalid_argument("DG::Integrator::step: dt must be > 0");

      // p_1 = u
      Real* du_data = du_.data().getData();
      Real* u_data = u.data().getData();
      Real* tmp_data = tmp_.data().getData();


      for (int i = 0; i < 5; i++)
      {
        // tmp = L(u, t_n + c_i dt)
        rhs_(u, tmp_, t_in + c_[i] * dt);

        // tmp = dt * tmp
        this -> scaleField_(tmp_, dt);

        // du_ = a_i * du_ + tmp
        this -> addScaled_(tmp_, du_, a_[i], du_);

        // u = u + b_i * du_
        this -> addScaled_(u, du_, b_[i], u);
      }
    }


    // integrate: perform multiple steps
    Index integrate(Field& u, Real t0, Real t_end, Real dt) override
    {
      if (t_end <= t0)
        throw std::invalid_argument("Integrator::integrate: t_end must be > t0");

      Real t = t0;
      Index nstep = 0;

      while (t < t_end - Real(1e-12))
      {
        // clip the last step to reach t_end exactly
        Real dt_actual = std::min(dt, t_end - t);

        step(u, dt_actual, t);
        t += dt_actual;
        nstep += 1;

        // Invoke optional callback; stop early if it returns false
        if (callback_ && !callback_(t, u, nstep))
            break;
      }

      return nstep;
    }

    // max dt from the Hesthaven book
    static Real computeDt(Real x_min, Real max_wave_speed, int poly_order, Real cfl = Real(0.4))
    {
        if (TNL::abs(max_wave_speed) <= Real(0))
            throw std::invalid_argument("max_wave_speed must be positive");
        return cfl * x_min / (TNL::abs(max_wave_speed) * ( 2 * poly_order + 1 ));
    }
    // getters
    Index numSteps() const override { return lastStepCount_; }
    Real currentTime() const override { return currentTime_; }
    Index numPoints() const override { return K_ * Np_; }

private:

    // RHS operator (function)
    RHSFunc rhs_;

    // number of cells and nodes
    Index K_, Np_;

    // the four stages
    Field du_;
    Field tmp_;

    // callback function
    Callback callback_;

    // coefficients, Hesthaven p. 64, or the NASA paper
    std::array<Real, 5> a_ = {
      0.f,
      -Real(567301805773.0)/Real(1357537059087.0),
      -Real(2404267990393.0)/Real(2016746695238.0),
      -Real(3550918686646.0)/Real(2091501179385.0),
      -Real(1275806237668.0)/Real(842570457699.0)
    };
    std::array<Real, 5> b_ = {
      Real(1432997174477.0)/Real(9575080441755.0),
      Real(5161836677717.0)/Real(13612068292357.0),
      Real(1720146321549.0)/Real(2090206949498.0),
      Real(3134564353537.0)/Real(4481467310338.0),
      Real(2277821191437.0)/Real(14882151754819.0)
    };
    std::array<Real, 5> c_ = {
      Real(0.f),
      Real(1432997174477.0)/Real(9575080441755.0),
      Real(2526269341429.0)/Real(6820363962896.0),
      Real(2006345519317.0)/Real(3224310063776.0),
      Real(2802321613138.0)/Real(2924317926251.0)
    };

    // diagnostics variables, updated by step()
    Index lastStepCount_{ 0 };
    Real currentTime_{ 0 };
};

// ------------------------------- SSP-RK3 ---------------------------------------
// strong stability preserving Runge-Kutta, four stage
template<typename Real   = double,
         typename Device = TNL::Devices::Host,
         typename Index  = int>
class SSPRK : Integrator< Real, Device, Index > // inherit protected methods
{
public:
    using Field = FieldVector<Real, Device, Index>;
    using RHSFunc = std::function<void(const Field&, Field&, const Real&)>;
    using Callback = StepCallback<Real, Device, Index>;

    // op       – spatial residual (must outlive this integrator)
    // K, Np    – element count and nodes-per-element (needed to size scratch)
    // callback – optional: called after every step, return false to stop early
    explicit SSPRK(RHSFunc rhs,
                 Index K, Index Np,
                 Callback callback = nullptr)
        : rhs_(std::move(rhs)),
          K_(K), Np_(Np),
          callback_(std::move(callback)),
          k1_(K, Np), k2_(K, Np),
          tmp_(K, Np) {}

    // single step
    void step(Field& u, Real dt, Real t_in) override
    {
      if (dt <= Real(0))
        throw std::invalid_argument("DG::Integrator::step: dt must be > 0");

      // k1 = u^n + 1/2 dt L(u^n, t_n)
      // tmp = L(u^n, t_n)
      rhs_(u, tmp_, t_in);
      // tmp = 1/2 dt tmp
      this -> scaleField_(tmp_, 0.5 * dt);
      // k1 = tmp + u^n
      this -> addScaled_(tmp_, u, 1.0, k1_);

      // k2 = k1 + 1/2 dt L(k1, t_n + 1/2 dt)
      // tmp = L(k1, t_n + 1/2 dt)
      rhs_(k1_, tmp_, t_in + 0.5 * dt);
      // tmp = 1/2 tmp
      this -> scaleField_(tmp_, 0.5 * dt);
      // k2 = tmp + u^n
      this -> addScaled_(tmp_, k1_, 1.0, k2_);

      // k3 = 2/3 u^n + 1/2 k2 + 1/6 L(k2, t^n + dt)
      // (k3 = k1) to save storage
      // tmp = L(k2, t^n + dt)
      rhs_(k2_, tmp_, t_in + dt);
      // tmp = 1/6 dt tmp
      this -> scaleField_(tmp_, Real(1.0/6.0) * dt);
      // tmp = tmp + 1/3 k2
      this -> addScaled_(tmp_, k2_, Real(1.0/3.0), tmp_);
      // k3 = tmp + 2/3 u^n
      this -> addScaled_(tmp_, u, Real(2.0/3.0), k1_);

      // u^n+1 = k3 + 1/2 dt L(k3, t^n + 1/2 dt)
      // tmp = L(k3, t^n + 1/2 dt)
      rhs_(k1_, tmp_, t_in + 0.5 * dt);
      // tmp = 1/2 dt tmp
      this -> scaleField_(tmp_, 0.5 * dt);
      // u = tmp + k3
      this -> addScaled_(tmp_, k1_, 1.0, u);

    }

    // integrate: perform multiple steps
    Index integrate(Field& u, Real t0, Real t_end, Real dt) override
    {
      if (t_end <= t0)
        throw std::invalid_argument("Integrator::integrate: t_end must be > t0");

      Real t = t0;
      Index nstep = 0;

      while (t < t_end - Real(1e-12))
      {
        // clip the last step to reach t_end exactly
        Real dt_actual = std::min(dt, t_end - t);

        step(u, dt_actual, t);
        t += dt_actual;
        nstep += 1;

        // Invoke optional callback; stop early if it returns false
        if (callback_ && !callback_(t, u, nstep))
            break;
      }

      return nstep;
    }

    // max dt from the Hesthaven book
    static Real computeDt(Real x_min, Real max_wave_speed, int poly_order, Real cfl = Real(0.4))
    {
        if (TNL::abs(max_wave_speed) <= Real(0))
            throw std::invalid_argument("max_wave_speed must be positive");
        return cfl * x_min / (TNL::abs(max_wave_speed) * ( 2 * poly_order + 1 ));
    }

    // getters
    Index numSteps() const override { return lastStepCount_; }
    Real currentTime() const override { return currentTime_; }
    Index numPoints() const override { return K_ * Np_; }

private:

    // RHS operator (function)
    RHSFunc rhs_;

    // number of cells and nodes
    Index K_, Np_;

    // the four stages
    Field k1_, k2_;
    Field tmp_;

    // callback function
    Callback callback_;

    // diagnostics variables, updated by step()
    Index lastStepCount_{ 0 };
    Real currentTime_{ 0 };
};

} // namespace DG
