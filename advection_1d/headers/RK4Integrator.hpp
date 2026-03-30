#pragma once

#include "FieldVector.hpp"
#include "Operator.hpp"
#include <functional>
#include <stdexcept>
#include <cmath>

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
class Integrator
{
public:
    using Field = FieldVector<Real, Device, Index>;
    using OperatorType = Operator<Real, Device, Index>;
    using Callback = StepCallback<Real, Device, Index>;

    virtual ~Integrator() = default;

    // single step
    virtual void step(Field& u, Real dt, Real t_in) = 0;

    // integrate, multiple steps
    virtual Index integrate(Field& u, Real t0, Real t_end, Real dt) = 0;

    // getters
    virtual Index numSteps() const = 0;
    virtual Index numPoints() const = 0;
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


// ------------------------------- RKSS ---------------------------------------
template<typename Real   = double,
         typename Device = TNL::Devices::Host,
         typename Index  = int>
class RKSS : Integrator< Real, Device, Index > // inherit protected methods
{
public:
    using Field = FieldVector<Real, Device, Index>;
    using OperatorType = Operator<Real, Device, Index>;
    using Callback = StepCallback<Real, Device, Index>;

    // op       – spatial residual (must outlive this integrator)
    // K, Np    – element count and nodes-per-element (needed to size scratch)
    // callback – optional: called after every step, return false to stop early
    explicit RKSS(const OperatorType& op,
                           Index K, Index Np,
                           Callback callback = nullptr)
        : op_(op)
        , K_(K), Np_(Np)
        , callback_(std::move(callback))
        , k1_(K, Np), k2_(K, Np), k3_(K, Np), k4_(K, Np)
        , tmp_(K, Np) {}

    // single step
    void step(Field& u, Real dt, Real t_in) override
    {
      if (dt <= Real(0))
        throw std::invalid_argument("DG::Integrator::step: dt must be > 0");

      // k1 = L(u^n, t_n)
      op_.computeRHS(u, k1_, t_in + 0.5 * dt);

      // k2 = L(u^n + 1/2 dt k1, t_n + 1/2 dt)
      // tmp = u + 0.5*dt*k1
      this -> addScaled_(u, k1_, 0.5 * dt, tmp_);
      op_.computeRHS(tmp_, k2_, t_in + 0.5 * dt);

      // k3 = L(u^n + 1/2 dt k2, t_n + 1/2 dt)
      // tmp = u + 0.5*dt*k2
      this -> addScaled_(u, k2_, 0.5 * dt, tmp_);
      op_.computeRHS(tmp_, k3_, t_in + 0.5 * dt);

      // k4 = L(u^n + dt k3, t_n + dt)
      // tmp = u + dt * k3
      this -> addScaled_(u, k3_, dt, tmp_);
      op_.computeRHS(tmp_, k4_, t_in + dt );

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
    static Real computeDt(Real x_min, Real max_wave_speed, Real cfl = Real(0.4))
    {
        if (max_wave_speed <= Real(0))
            throw std::invalid_argument("max_wave_speed must be positive");
        return cfl * x_min / max_wave_speed;
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

    // RHS operator
    const OperatorType& op_;

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

} // namespace DG
