// RK4Integrator.hpp
#pragma once
#include "FieldVector.hpp"
#include "Operator.hpp"
#include <functional>
#include <stdexcept>
#include <cmath>

namespace DG {

// ----------------------------------------------------------------------------
//  Callback type: called after every accepted step with (t, u, step_index).
//  Use it for output, diagnostics, or early stopping (return false to halt).
// ----------------------------------------------------------------------------
template<typename Real, typename Device, typename Index>
using StepCallback = std::function< bool(Real t, const FieldVector<Real,Device,Index>& u, Index step) >;

// ----------------------------------------------------------------------------
//  RK4Integrator
//
//  Advances a DGFieldVector in time using the classical 4-stage, 4th-order
//  Runge-Kutta method, i.e. the "RK4" formula from every ODE textbook:
//
//      k1 = dt * L(u^n)
//      k2 = dt * L(u^n + k1/2)
//      k3 = dt * L(u^n + k2/2)
//      k4 = dt * L(u^n + k3)
//      u^{n+1} = u^n + (k1 + 2*k2 + 2*k3 + k4) / 6
//
//  where L = DGOperator::computeRHS is the spatial residual.
//
//  Template parameters mirror DGFieldVector and DGOperator:
//    Real   – floating-point type (double or float)
//    Device – TNL device (TNL::Devices::Host or TNL::Devices::Cuda)
//    Index  – integer type for element/DOF counts
// ----------------------------------------------------------------------------
template<typename Real   = double,
         typename Device = TNL::Devices::Host,
         typename Index  = int>
class Integrator
{
public:
    using Field    = FieldVector<Real, Device, Index>;
    using Operator = Operator<Real, Device, Index>;
    using Callback = StepCallback<Real, Device, Index>;

    // ------------------------------------------------------------------ //
    //  Construction
    // ------------------------------------------------------------------ //

    // op       – spatial residual (must outlive this integrator)
    // K, Np    – element count and nodes-per-element (needed to size scratch)
    // callback – optional: called after every step, return false to stop early
    explicit Integrator(const Operator& op,
                           Index K, Index Np,
                           Callback callback = nullptr)
        : op_(op)
        , K_(K), Np_(Np)
        , callback_(std::move(callback))
        , k1_(K, Np), k2_(K, Np), k3_(K, Np), k4_(K, Np)
        , tmp_(K, Np)
    {}

    // ------------------------------------------------------------------ //
    //  Single step:  u ← u + RK4(u, dt)
    //
    //  Advances the solution in-place by one time step of size dt.
    //  t_in is the time at the start of this step (only used for the
    //  callback; the operator itself is autonomous and does not need t).
    // ------------------------------------------------------------------ //
    void step(Field& u, Real dt, Real t_in = Real(0))
    {
        if (dt <= Real(0))
            throw std::invalid_argument("RK4Integrator::step: dt must be > 0");

        // --- stage 1: k1 = dt * L(u^n) ---
        op_.computeRHS(u, k1_);
        scaleField_(k1_, dt);                       // k1 *= dt

        // --- stage 2: k2 = dt * L(u^n + k1/2) ---
        addScaled_(u, k1_, Real(0.5), tmp_);        // tmp = u + 0.5*k1
        op_.computeRHS(tmp_, k2_);
        scaleField_(k2_, dt);                       // k2 *= dt

        // --- stage 3: k3 = dt * L(u^n + k2/2) ---
        addScaled_(u, k2_, Real(0.5), tmp_);        // tmp = u + 0.5*k2
        op_.computeRHS(tmp_, k3_);
        scaleField_(k3_, dt);                       // k3 *= dt

        // --- stage 4: k4 = dt * L(u^n + k3) ---
        addScaled_(u, k3_, Real(1.0), tmp_);        // tmp = u + k3
        op_.computeRHS(tmp_, k4_);
        scaleField_(k4_, dt);                       // k4 *= dt

        // --- combine: u += (k1 + 2*k2 + 2*k3 + k4) / 6 ---
        combine_(u, k1_, k2_, k3_, k4_);
    }

    // ------------------------------------------------------------------ //
    //  Time integration loop: integrate from t0 to t_end with step dt.
    //
    //  The last step is clipped to land exactly on t_end (avoids
    //  overshooting due to floating-point remainder).
    //
    //  Returns the number of steps taken.
    // ------------------------------------------------------------------ //
    Index integrate(Field& u, Real t0, Real t_end, Real dt)
    {
        if (t_end <= t0)
            throw std::invalid_argument("RK4Integrator::integrate: t_end must be > t0");
        if (dt <= Real(0))
            throw std::invalid_argument("RK4Integrator::integrate: dt must be > 0");

        Real t      = t0;
        Index nstep = 0;

        while (t < t_end - Real(1e-12) * std::abs(t_end))
        {
            // Clip the last step to reach t_end exactly
            Real dt_actual = std::min(dt, t_end - t);

            step(u, dt_actual, t);
            t     += dt_actual;
            nstep += 1;

            // Invoke optional callback; stop early if it returns false
            if (callback_ && !callback_(t, u, nstep))
                break;
        }

        return nstep;
    }

    // ------------------------------------------------------------------ //
    //  Convenience: compute a stable dt given mesh and polynomial order.
    //
    //  Uses the standard DG CFL estimate:
    //      dt = CFL * h_min / (|a| * (2N+1)^2)
    //
    //  cfl should be in (0, 1]; 0.4 is a safe default for RK4 + DG.
    // ------------------------------------------------------------------ //
    static Real computeDt(Real h_min, Real max_wave_speed,
                          Index poly_order, Real cfl = Real(0.4))
    {
        if (max_wave_speed <= Real(0))
            throw std::invalid_argument("max_wave_speed must be positive");
        return cfl * h_min / (max_wave_speed
                               * std::pow(Real(2 * poly_order + 1), Real(2)));
    }

    // ------------------------------------------------------------------ //
    //  Accessors
    // ------------------------------------------------------------------ //
    Index numSteps()    const { return lastStepCount_; }
    Real  currentTime() const { return currentTime_;   }

private:
    // ------------------------------------------------------------------ //
    //  Helper: scale every DOF of f by scalar s   (f *= s)
    // ------------------------------------------------------------------ //
    void scaleField_(Field& f, Real s) const
    {
        const Index total = K_ * Np_;
        Real* data = f.data().getData();
        for (Index i = 0; i < total; ++i) data[i] *= s;
    }

    // ------------------------------------------------------------------ //
    //  Helper: out = base + alpha * delta   (element-wise)
    // ------------------------------------------------------------------ //
    void addScaled_(const Field& base, const Field& delta,
                    Real alpha, Field& out) const
    {
        const Index total = K_ * Np_;
        const Real* b = base.data().getData();
        const Real* d = delta.data().getData();
        Real*       o = out.data().getData();
        for (Index i = 0; i < total; ++i) o[i] = b[i] + alpha * d[i];
    }

    // ------------------------------------------------------------------ //
    //  Helper: u += (k1 + 2*k2 + 2*k3 + k4) / 6
    // ------------------------------------------------------------------ //
    void combine_(Field& u,
                  const Field& k1, const Field& k2,
                  const Field& k3, const Field& k4) const
    {
        const Index total = K_ * Np_;
        Real*       pu = u.data().getData();
        const Real* p1 = k1.data().getData();
        const Real* p2 = k2.data().getData();
        const Real* p3 = k3.data().getData();
        const Real* p4 = k4.data().getData();
        const Real sixth = Real(1) / Real(6);
        for (Index i = 0; i < total; ++i)
            pu[i] += sixth * (p1[i] + Real(2)*p2[i] + Real(2)*p3[i] + p4[i]);
    }

    // ------------------------------------------------------------------ //
    //  Data members
    // ------------------------------------------------------------------ //
    const Operator& op_;           // reference to the spatial residual
    Index           K_, Np_;       // mesh and basis sizes
    Callback        callback_;     // optional post-step hook

    // Scratch fields: allocated once in constructor, reused every step
    Field k1_, k2_, k3_, k4_;     // the four RK4 stage derivatives
    Field tmp_;                    // temporary u + alpha*k_i

    // Diagnostics (updated by integrate())
    Index lastStepCount_{ 0 };
    Real  currentTime_  { 0 };
};

} // namespace DG
