#include <gtest/gtest.h>
#include "../headers/Mesh.hpp"
#include "../headers/ReferenceElement.hpp"
#include "../headers/FieldVector.hpp"
#include "../headers/NumericalFlux.hpp"
#include "../headers/Operator.hpp"
#include "../headers/RK4Integrator.hpp"
#include <cmath>

TEST(RK4Test, AdvectsWaveOneFullPeriod)
{
    // After T = 2π/a = 2π steps, sin(x - at) returns to sin(x).
    // The L2 error should be small for N=4, K=10.
    const int K = 10, N = 4;
    const double a = 1.0;
    const double Tf = 2.0 * M_PI;

    auto mesh = DG::Mesh<double>::uniform(0.0, Tf, K);
    DG::ReferenceElement<double> ref(N);
    DG::UpwindFlux<double>       flux(a);
    DG::Operator<double>         op(mesh, ref, flux,
                                    [](double u){ return u; }, a);

    int Np = ref.numDOF();
    DG::FieldVector<double> u(K, Np);

    // u(x,0) = sin(x)
    for (int k = 0; k < K; ++k) {
        double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
        for (int i = 0; i < Np; ++i) {
            double r = ref.nodes()[i];
            u.elementPtr(k)[i] = std::sin(xL + (r+1)*0.5*h);
        }
    }

    double h_min = mesh.elementSize(0);
    double dt    = DG::RK4Integrator<double>::computeDt(h_min, a, N);

    DG::RK4Integrator<double> rk(op, K, Np);
    rk.integrate(u, 0.0, Tf, dt);

    // Compute L2 error against sin(x) (exact solution at T=2π)
    double err = 0;
    for (int k = 0; k < K; ++k) {
        double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
        for (int i = 0; i < Np; ++i) {
            double r = ref.nodes()[i];
            double x = xL + (r+1)*0.5*h;
            double e = u.elementPtr(k)[i] - std::sin(x);
            err += ref.weights()[i] * mesh.jacobian(k) * e * e;
        }
    }
    err = std::sqrt(err);

    // N=4 DG should achieve L2 error well below 1e-3 for K=10
    EXPECT_LT(err, 1e-3) << "L2 error after one period = " << err;
}

TEST(RK4Test, EnergyConservedOverTime)
{
    const int K = 10, N = 4;
    const double a = 1.0;
    const double Tf = 1.0;

    auto mesh = DG::Mesh<double>::uniform(0.0, 2*M_PI, K);
    DG::ReferenceElement<double> ref(N);
    DG::UpwindFlux<double>       flux(a);
    DG::Operator<double>         op(mesh, ref, flux,
                                    [](double u){ return u; }, a);

    int Np = ref.numDOF();
    DG::FieldVector<double> u(K, Np);

    for (int k = 0; k < K; ++k) {
        double xL = mesh.leftVertex(k), h = mesh.elementSize(k);
        for (int i = 0; i < Np; ++i) {
            double r = ref.nodes()[i];
            u.elementPtr(k)[i] = std::sin(xL + (r+1)*0.5*h);
        }
    }

    auto energy = [&](const DG::FieldVector<double>& v) {
        double E = 0;
        for (int k = 0; k < K; ++k) {
            double J = mesh.jacobian(k);
            for (int i = 0; i < Np; ++i)
                E += ref.weights()[i] * J * v.elementPtr(k)[i] * v.elementPtr(k)[i];
        }
        return E;
    };

    double E0 = energy(u);
    double h_min = mesh.elementSize(0);
    double dt    = DG::RK4Integrator<double>::computeDt(h_min, a, N);

    DG::RK4Integrator<double> rk(op, K, Np);
    rk.integrate(u, 0.0, Tf, dt);

    double E1 = energy(u);
    // Energy should not grow — allow small decrease from upwind dissipation
    EXPECT_LE(E1, E0 * 1.001) << "Energy grew from " << E0 << " to " << E1;
}
