#pragma once

#include "FieldVector.hpp"
#include <algorithm>
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
class Integrator
{
public:
    using Field = FieldVector<Real, Device, Index>;
    using RHSFunc = std::function<void(const Field&, Field&, const Real&)>;
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


// ------------------------------- ERK ---------------------------------------
// simple explicit four stage Runge-Kutta
template<typename Real   = double,
         typename Device = TNL::Devices::Host,
         typename Index  = int>
class ERK : Integrator< Real, Device, Index > // inherit protected methods
{
public:
    using Field = FieldVector<Real, Device, Index>;
    using RHSFunc = std::function<void(const Field&, Field&, const Real&)>;
    using Callback = StepCallback<Real, Device, Index>;

    // op       – spatial residual (must outlive this integrator)
    // K, Np    – element count and nodes-per-element (needed to size scratch)
    // callback – optional: called after every step, return false to stop early
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
        , k1_(K, Np), k2_(K, Np), p1_(K, Np), p2_(K, Np)
        , tmp_(K, Np) {}

    // single step
    void step(Field& u, Real dt, Real t_in) override
    {
      if (dt <= Real(0))
        throw std::invalid_argument("DG::Integrator::step: dt must be > 0");

      // p_1 = u
      for (int i = 0; i < K_ * Np_; i++)
        p1_[i] = u[i];

      for (int i = 0; i < 5; i++)
      {
        // tmp = L(p1, t_n + c_i dt)
        rhs_(p1_, tmp_, t_in + c_[i] * dt);
        // k_i = a_i k_(i-1) + tmp
        this -> addScaled_(tmp_, k1_, a_[i], k2_);
        // p_i = p_(i-1) + b_i * k_i
        this -> addScaled_(p1_, k2_, b_[i], p2_);
      }

      // u = p_5
      for (int i = 0; i < K_ * Np_; i++)
        u[i] = p2_[i];
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
    
    // RHS operator (function)
    RHSFunc rhs_;

    // number of cells and nodes
    Index K_, Np_;

    // the four stages
    Field k1_, k2_, p1_, p2_;
    Field tmp_;

    // callback function
    Callback callback_;

    // coefficients, Hesthaven p. 64
    std::vector<Real> a_ = {
      0.f,
      -567301805773/1357537059087,
      -2404267990393/2016746695238,
      -3550918686646/2091501179385,
      -1275806237668/842570457699
    };
    std::vector<Real> b_ = {
      1432997174477/9575080441755,
      5161836677717/13612068292357,
      1720146321549/2090206949498,
      3134564353537/4481467310338,
      2277821191437/14882151754819
    };
    std::vector<Real> c_ = {
      0.f,
      1432997174477/9575080441755,
      2526269341429/6820363962896,
      2006345519317/3224310063776,
      2802321613138/2924317926251
    };
    
    // diagnostics variables, updated by step()
    Index lastStepCount_{ 0 };
    Real currentTime_{ 0 };
};

} // namespace DG
