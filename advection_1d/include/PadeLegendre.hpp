/**
 * @file PadeLegendre.hpp
 * @brief Padé–Legendre post-processing reconstructor for DG solutions.
 *
 * Given element-local nodal DG values, builds an [L/M] Padé–Legendre
 * rational approximant P_L(x)/Q_M(x) on the reference element [-1,1].
 * This replaces the local polynomial with a rational function that can
 * better resolve sharp features and near-discontinuities.
 *
 * Mathematical background
 * -----------------------
 * On each element the DG solution is expressed as
 *   u_h(x) = sum_{n=0}^{N} c_n P_n(x)
 * where c_n are the (unnormalized) Legendre modal coefficients obtained via
 * the inverse Vandermonde: c = diag(sqrt((2n+1)/2)) * V^{-1} * u_nodal.
 *
 * The [L/M] Padé–Legendre approximant seeks
 *   P_L(x) = sum_{l=0}^{L} p_l P_l(x)   (numerator, degree L)
 *   Q_M(x) = 1 + sum_{m=1}^{M} q_m P_m(x)  (denominator, degree M, q_0=1)
 * such that the first N+1 Legendre moments of Q_M * u_h - P_L vanish.
 *
 * This yields:
 *   - An (N-L) × M denominator system  A q = -b  solved least-squares
 *   - An explicit formula  p_k = (2k+1)/2 * sum_m q_m B(k,m)  for k=0..L
 * where B(k,m) = integral_{-1}^{1} P_k(x) P_m(x) u_h(x) dx is computed
 * numerically with a Gauss–Legendre quadrature of order N+M+1.
 */

#pragma once

#include "FieldVector.hpp"
#include "ReferenceElement.hpp"
#include <Eigen/Dense>
#include <TNL/Devices/Host.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace DG {

// ─────────────────────────────────────────────────────────────────────────────
//  PadeApproximant
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Stores the [L/M] Padé–Legendre approximant for a single element.
 *
 * Coefficients are in the unnormalized Legendre polynomial basis.
 * q_coeffs[0] is always 1 (normalization convention).
 */
template<class Real = double, class Index = int>
struct PadeApproximant
{
    std::vector<Real> p_coeffs;     ///< Numerator  P_L: p[l] for l = 0..L
    std::vector<Real> q_coeffs;     ///< Denominator Q_M: q[m] for m = 0..M (q[0]==1)
    std::vector<Real> modal_coeffs; ///< u_h modal coefficients c_n (polynomial fallback)
    bool  valid{true};              ///< False if the linear solve failed
    Index L{0};                     ///< Numerator degree
    Index M{0};                     ///< Denominator degree
};

// ─────────────────────────────────────────────────────────────────────────────
//  PadeLegendreSolver
// ─────────────────────────────────────────────────────────────────────────────

/**
 * @class PadeLegendreSolver
 * @brief Builds and evaluates element-local [L/M] Padé–Legendre approximants.
 *
 * Constructed once per polynomial order; then used to process any number of
 * element-local nodal value arrays.  All arithmetic is done on the reference
 * element [-1,1] using unnormalized Legendre polynomials.
 */
template<class Real = double, class Index = int>
class PadeLegendreSolver
{
public:
    using RefElem   = ReferenceElement<Real, Index>;
    using RefMatrix = typename RefElem::Matrix;
    using Field     = FieldVector<Real, TNL::Devices::Host, Index>;
    using EigenMat  = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using EigenVec  = Eigen::Matrix<double, Eigen::Dynamic, 1>;

    /**
     * @brief Construct the solver.
     *
     * Precomputes the inverse Vandermonde matrix and a Gauss–Legendre
     * quadrature rule accurate enough to integrate P_k P_m u_h exactly.
     *
     * @param ref          Reference element (must outlive this object)
     * @param L            Numerator degree  (L >= 0)
     * @param M            Denominator degree (M >= 0, L+M <= N)
     * @param fallback_tol Absolute tolerance: if |Q(x)| < tol use polynomial fallback.
     *                     Defaults to 100 * machine epsilon for type Real, so
     *                     the value is automatically appropriate for both float
     *                     and double instantiations.
     */
    PadeLegendreSolver(const RefElem& ref, Index L, Index M,
                       Real fallback_tol = Real(100) * std::numeric_limits<Real>::epsilon())
        : ref_(ref), L_(L), M_(M), N_(ref.order()), Np_(ref.numDOF()),
          fallback_tol_(fallback_tol), vinv_(ref.Vinv())
    {
        if (L < 0 || M < 0)
            throw std::invalid_argument("PadeLegendreSolver: L and M must be >= 0");
        if (L + M > N_)
            throw std::invalid_argument("PadeLegendreSolver: L+M must be <= N");

        // Need 2*nq-1 >= k_max + m_max + N = N + M + N = 2N + M
        // nq = N + M + 1 gives 2(N+M+1)-1 = 2N+2M+1 >= 2N+M  ✓
        const Index nq = N_ + M_ + 1;
        buildGaussLegendre_(nq, gl_nodes_, gl_weights_);
    }

    // ── Public interface ──────────────────────────────────────────────────────

    /**
     * @brief Convert element-local nodal values to unnormalized modal coefficients.
     *
     * Returns c such that u_h(x) = sum_{n=0}^{N} c[n] * P_n(x).
     * The Vandermonde inverse gives normalized coefficients c_hat_n; these are
     * converted by c_n = c_hat_n * sqrt((2n+1)/2).
     *
     * @param nodal_vals Pointer to Np nodal values at the GLL nodes
     * @return vector of Np unnormalized Legendre coefficients
     */
    std::vector<Real> modalCoeffs(const Real* nodal_vals) const
    {
        std::vector<Real> c(Np_);
        for (Index n = 0; n < Np_; ++n) {
            Real s = Real(0);
            for (Index i = 0; i < Np_; ++i)
                s += vinv_(n, i) * nodal_vals[i];
            // vinv_ applied to nodal gives normalized coefficient c_hat_n = s
            // unnormalize: c_n = c_hat_n * sqrt((2n+1)/2)
            c[n] = s * std::sqrt(Real(2 * n + 1) / Real(2));
        }
        return c;
    }

    /**
     * @brief Build the [L/M] Padé–Legendre approximant from modal coefficients.
     *
     * Solves the denominator system A q = -b (least-squares via
     * colPivHouseholderQr), then computes numerator coefficients explicitly.
     *
     * @param c Modal coefficients from modalCoeffs()
     * @return PadeApproximant containing p_coeffs, q_coeffs and modal_coeffs
     */
    PadeApproximant<Real, Index> buildApproximant(const std::vector<Real>& c) const
    {
        PadeApproximant<Real, Index> approx;
        approx.L = L_;
        approx.M = M_;
        approx.modal_coeffs = c;
        approx.p_coeffs.assign(L_ + 1, Real(0));
        approx.q_coeffs.assign(M_ + 1, Real(0));
        approx.q_coeffs[0] = Real(1);  // normalization

        // B(k,m) = integral_{-1}^1 P_k(x) P_m(x) u_h(x) dx
        const EigenMat B = buildBMatrix_(c);

        if (M_ == 0) {
            // Denominator is trivially 1; numerator equals polynomial truncation
            for (Index k = 0; k <= L_; ++k)
                approx.p_coeffs[k] = static_cast<Real>(
                    double(2 * k + 1) / 2.0 * B(k, 0));
            approx.valid = true;
            return approx;
        }

        // ── Denominator system: A q_vec = -b ────────────────────────────────
        // Equations for k = L+1 .. N  →  (N-L) equations, M unknowns
        const Index nRows = N_ - L_;
        EigenMat A(nRows, M_);
        EigenVec b(nRows);

        for (Index k = L_ + 1; k <= N_; ++k) {
            const Index row = k - L_ - 1;
            b(row) = B(k, 0);
            for (Index m = 1; m <= M_; ++m)
                A(row, m - 1) = B(k, m);
        }

        // Least-squares solve (handles square and overdetermined systems)
        const EigenVec q_vec = A.colPivHouseholderQr().solve(-b);

        for (Index m = 1; m <= M_; ++m)
            approx.q_coeffs[m] = static_cast<Real>(q_vec(m - 1));

        // ── Numerator: p_k = (2k+1)/2 * sum_{m=0}^M q_m B(k,m) ─────────────
        for (Index k = 0; k <= L_; ++k) {
            double s = 0.0;
            for (Index m = 0; m <= M_; ++m)
                s += static_cast<double>(approx.q_coeffs[m]) * B(k, m);
            approx.p_coeffs[k] = static_cast<Real>(double(2 * k + 1) / 2.0 * s);
        }

        approx.valid = true;
        return approx;
    }

    /**
     * @brief Evaluate the Padé approximant P(x)/Q(x) at reference coordinate x.
     *
     * If |Q(x)| < fallback_tol, returns the polynomial approximation u_h(x)
     * instead to avoid division by a near-zero denominator.
     *
     * @param approx Previously built PadeApproximant
     * @param x      Reference coordinate in [-1,1]
     * @return Reconstructed value at x
     */
    Real evaluate(const PadeApproximant<Real, Index>& approx, Real x) const
    {
        Real P = Real(0);
        for (Index l = 0; l <= approx.L; ++l)
            P += approx.p_coeffs[l] * RefElem::legendreP(l, x);

        Real Q = Real(0);
        for (Index m = 0; m <= approx.M; ++m)
            Q += approx.q_coeffs[m] * RefElem::legendreP(m, x);

        if (std::abs(Q) < fallback_tol_) {
            // Polynomial fallback: evaluate the original DG approximation
            Real poly = Real(0);
            for (Index n = 0; n < static_cast<Index>(approx.modal_coeffs.size()); ++n)
                poly += approx.modal_coeffs[n] * RefElem::legendreP(n, x);
            return poly;
        }

        return P / Q;
    }

    /**
     * @brief Reconstruct at a vector of reference coordinates.
     *
     * Convenience wrapper: converts nodal values → modal coefficients →
     * approximant → evaluations at output_nodes.
     *
     * @param nodal_vals  Pointer to Np nodal values
     * @param output_nodes Reference coordinates at which to evaluate
     * @return Reconstructed values at each output node
     */
    std::vector<Real> reconstruct(const Real*               nodal_vals,
                                   const std::vector<Real>&  output_nodes) const
    {
        const auto c      = modalCoeffs(nodal_vals);
        const auto approx = buildApproximant(c);
        std::vector<Real> out(output_nodes.size());
        for (std::size_t i = 0; i < output_nodes.size(); ++i)
            out[i] = evaluate(approx, output_nodes[i]);
        return out;
    }

    /**
     * @brief Apply element-wise reconstruction to a full FieldVector.
     *
     * Each element is treated independently on the reference element [-1,1].
     * Returns a new FieldVector with the reconstructed values at the same GLL
     * nodes as the input.
     *
     * @param u Input DG solution FieldVector
     * @return  New FieldVector with Padé-reconstructed nodal values
     */
    Field reconstruct(const Field& u) const
    {
        Field result(u.numElements(), Np_);

        // Collect GLL node positions once
        std::vector<Real> gll_nodes(Np_);
        for (Index i = 0; i < Np_; ++i)
            gll_nodes[i] = ref_.nodes()[i];

        for (Index k = 0; k < u.numElements(); ++k) {
            const Real* uk = u.elementPtr(k);
            const auto  vals = reconstruct(uk, gll_nodes);
            Real* rk = result.elementPtr(k);
            for (Index i = 0; i < Np_; ++i)
                rk[i] = vals[i];
        }
        return result;
    }

private:

    // ── Private helpers ───────────────────────────────────────────────────────

    /**
     * @brief Build the B matrix using Gauss–Legendre quadrature.
     *
     * B(k,m) = integral_{-1}^{1} P_k(x) P_m(x) u_h(x) dx
     * for k = 0..N, m = 0..M.
     *
     * @param c Unnormalized modal coefficients of u_h
     * @return  (N+1) x (M+1) Eigen matrix
     */
    EigenMat buildBMatrix_(const std::vector<Real>& c) const
    {
        const Index nq = static_cast<Index>(gl_nodes_.size());
        EigenMat B = EigenMat::Zero(N_ + 1, M_ + 1);

        for (Index q = 0; q < nq; ++q) {
            const Real xq = gl_nodes_[q];
            const Real wq = gl_weights_[q];

            // Evaluate u_h at quadrature point
            Real uh = Real(0);
            for (Index n = 0; n <= N_; ++n)
                uh += c[n] * RefElem::legendreP(n, xq);

            // Accumulate B(k, m) += w_q * P_k(x_q) * P_m(x_q) * u_h(x_q)
            for (Index k = 0; k <= N_; ++k) {
                const double Pk = static_cast<double>(RefElem::legendreP(k, xq));
                for (Index m = 0; m <= M_; ++m) {
                    const double Pm = static_cast<double>(RefElem::legendreP(m, xq));
                    B(k, m) += static_cast<double>(wq) * Pk * Pm
                               * static_cast<double>(uh);
                }
            }
        }
        return B;
    }

    /**
     * @brief Compute Gauss–Legendre quadrature on [-1,1] via Golub–Welsch.
     *
     * Builds the symmetric tridiagonal Jacobi matrix for the Legendre
     * polynomial family, then finds its eigendecomposition.  Eigenvalues
     * are the quadrature nodes; weights are 2*(first eigenvector entry)^2.
     *
     * @param nq      Number of quadrature points
     * @param nodes   Output: quadrature nodes
     * @param weights Output: quadrature weights
     */
    void buildGaussLegendre_(Index nq,
                              std::vector<Real>& nodes,
                              std::vector<Real>& weights) const
    {
        if (nq < 1)
            throw std::invalid_argument("PadeLegendreSolver: nq must be >= 1");

        // Symmetric tridiagonal Jacobi matrix for Gauss–Legendre:
        //   off-diagonal entry: beta_k = k / sqrt(4k^2 - 1)
        EigenMat J = EigenMat::Zero(nq, nq);
        for (int k = 1; k < nq; ++k) {
            const double beta = static_cast<double>(k)
                                / std::sqrt(static_cast<double>(4 * k * k - 1));
            J(k - 1, k) = beta;
            J(k, k - 1) = beta;
        }

        Eigen::SelfAdjointEigenSolver<EigenMat> solver(J);
        const auto& eigenvals = solver.eigenvalues();
        const auto& eigenvecs = solver.eigenvectors();

        nodes.resize(nq);
        weights.resize(nq);
        for (int i = 0; i < nq; ++i) {
            nodes[i]   = static_cast<Real>(eigenvals(i));
            const double v0 = eigenvecs(0, i);
            weights[i] = static_cast<Real>(2.0 * v0 * v0);
        }
    }

    // ── Member data ──────────────────────────────────────────────────────────

    const RefElem& ref_;          ///< Reference element (borrowed)
    Index L_, M_, N_, Np_;        ///< Padé degrees and element polynomial order
    Real  fallback_tol_;          ///< Denominator zero-tolerance for fallback
    std::vector<Real> gl_nodes_;  ///< Gauss–Legendre quadrature nodes
    std::vector<Real> gl_weights_;///< Gauss–Legendre quadrature weights
    RefMatrix vinv_;              ///< Cached inverse Vandermonde (normalized basis)
};

} // namespace DG
