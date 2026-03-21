#include "IO.hpp"
#include "RK4Integrator.hpp"
#include <iostream>
#include <cmath>

int main(int argc, char* argv[])
{
    using Real = double;
    const int  K = 10, N = 4;
    const Real a = 1.0, Tf = 2.0 * M_PI;

    // Build mesh: either from a VTK file or uniform
    DG::Mesh<Real> mesh = (argc > 1)
        ? DG::Mesh<Real>::readVTK(argv[1])
        : DG::Mesh<Real>::uniform(0.0, 2.0 * M_PI, K);

    DG::ReferenceElement<Real> ref(N);
    DG::UpwindFlux<Real> flux(a);
    DG::Operator<Real> op(mesh, ref, flux, [&](Real u){ return a * u; });

    const int Np = ref.numDOF();
    DG::FieldVector<Real> u(mesh.numElements(), Np);

    // Initial condition u(x,0) = sin(x)
    for (int k = 0; k < mesh.numElements(); ++k) {
        Real xL = mesh.leftVertex(k), h = mesh.elementSize(k);
        Real* uk = u.elementPtr(k);
        for (int i = 0; i < Np; ++i) {
            Real r = ref.nodes()[i];
            uk[i]  = std::sin(xL + (r + 1.0) * 0.5 * h);
        }
    }

    // Write initial condition
    int frame = 0;
    DG::writeTimeSeriesVTK(mesh, ref, u, "output", frame++, Real(0));

    // Time stepping with output every 20 steps
    Real h_min = mesh.elementSize(0);
    Real dt    = DG::Integrator<Real>::computeDt(h_min, a, N);

//#if DIAGNOSTICS
    // --- diagnostics ---
    std::cout << "h_min = " << h_min << "  dt = " << dt << "\n\n";

    std::cout << "GLL nodes:   ";
    for (int i = 0; i < ref.numDOF(); ++i) std::cout << ref.nodes()[i] << " ";
    std::cout << "\n";

    std::cout << "GLL weights: ";
    for (int i = 0; i < ref.numDOF(); ++i) std::cout << ref.weights()[i] << " ";
    std::cout << "\n";

    // sum of weights should equal 2.0 (length of [-1,1])
    Real wsum = 0;
    for (int i = 0; i < ref.numDOF(); ++i) wsum += ref.weights()[i];
    std::cout << "sum(weights) = " << wsum << "  (should be 2.0)\n\n";

    // D matrix: D*ones should be zero (derivative of constant = 0)
    std::cout << "D * ones (should all be ~0):\n";
    for (int i = 0; i < ref.numDOF(); ++i) {
        Real s = 0;
        for (int j = 0; j < ref.numDOF(); ++j) s += ref.Dr()(i,j);
        std::cout << "  row " << i << ": " << s << "\n";
    }

    // D matrix: D*r should be ones (derivative of identity = 1)
    std::cout << "D * r (should all be ~1.0):\n";
    for (int i = 0; i < ref.numDOF(); ++i) {
        Real s = 0;
        for (int j = 0; j < ref.numDOF(); ++j) s += ref.Dr()(i,j) * ref.nodes()[j];
        std::cout << "  row " << i << ": " << s << "\n";
    }

    // LIFT: columns should sum to Np (surface-to-volume consistency)
    std::cout << "\nLIFT col 0 (left face):  ";
    Real lsum0 = 0;
    for (int i = 0; i < ref.numDOF(); ++i) { std::cout << ref.LIFT()(i,0) << " "; lsum0 += ref.LIFT()(i,0); }
    std::cout << "\n  sum = " << lsum0 << "\n";

    std::cout << "LIFT col 1 (right face): ";
    Real lsum1 = 0;
    for (int i = 0; i < ref.numDOF(); ++i) { std::cout << ref.LIFT()(i,1) << " "; lsum1 += ref.LIFT()(i,1); }
    std::cout << "\n  sum = " << lsum1 << "\n";
    // end of diagnostics
    //
    // 2. Compute RHS at t=0 and print element 0
    DG::FieldVector<Real> rhs0(mesh.numElements(), Np);
    op.computeRHS(u, rhs0);
    std::cout << "u[0]:   ";
    for (int i = 0; i < Np; ++i) std::cout << u.elementPtr(0)[i] << " ";
    std::cout << "\nrhs[0]: ";
    for (int i = 0; i < Np; ++i) std::cout << rhs0.elementPtr(0)[i] << " ";
    std::cout << "\n\n";

    // 3. For linear advection u_t + u_x = 0, exact rhs = -du/dx.
    // Compute -du/dx at element 0 manually using D matrix.
    std::cout << "exact rhs[0] = -J^{-1} * D * u[0]:\n";
    Real Jinv = 1.0 / mesh.jacobian(0);
    for (int i = 0; i < Np; ++i) {
        Real s = 0;
        for (int j = 0; j < Np; ++j)
            s += ref.Dr()(i,j) * u.elementPtr(0)[j];
        std::cout << "  row " << i << ": " << -Jinv * s << "\n";
    }
    // Print what the operator sees at each face of element 0
    // Face 0 (left face of element 0): left neighbour = last element (periodic)
    int kLast = mesh.numElements() - 1;
    std::cout << "Element 0 left face:\n";
    std::cout << "  u_minus (right end of last element): "
              << u.elementPtr(kLast)[Np-1] << "\n";
    std::cout << "  u_plus  (left end of element 0):     "
              << u.elementPtr(0)[0] << "\n";

    // Face 1 (right face of element 0): right neighbour = element 1
    std::cout << "Element 0 right face:\n";
    std::cout << "  u_minus (right end of element 0): "
              << u.elementPtr(0)[Np-1] << "\n";
    std::cout << "  u_plus  (left end of element 1):  "
              << u.elementPtr(1)[0] << "\n";

    // Also print u of element 1 left node and last element right node
    std::cout << "\nAll element left/right boundary values:\n";
    for (int k = 0; k < mesh.numElements(); ++k) {
        std::cout << "  k=" << k
                  << "  u_left="  << u.elementPtr(k)[0]
                  << "  u_right=" << u.elementPtr(k)[Np-1] << "\n";
    }

    // Manual forward Euler step and energy check
    op.computeRHS(u, rhs0);

    // Energy before: E = sum_k sum_i w_i * J_k * u_ki^2
    auto energy = [&](const DG::FieldVector<Real>& v) {
        Real E = 0;
        for (int k = 0; k < mesh.numElements(); ++k) {
            Real J = mesh.jacobian(k);
            const Real* vk = v.elementPtr(k);
            for (int i = 0; i < Np; ++i)
                E += ref.weights()[i] * J * vk[i] * vk[i];
        }
        return E;
    };

    std::cout << "Energy at t=0:         " << energy(u) << "\n";

    // Take one Euler step manually
    DG::FieldVector<Real> u1(mesh.numElements(), Np);
    for (int k = 0; k < mesh.numElements(); ++k)
        for (int i = 0; i < Np; ++i)
            u1.elementPtr(k)[i] = u.elementPtr(k)[i] + dt * rhs0.elementPtr(k)[i];

    std::cout << "Energy after 1 Euler step: " << energy(u1) << "\n";

    // Now compute rhs at u1 and check element 0 flux jump
    DG::FieldVector<Real> rhs1(mesh.numElements(), Np);
    op.computeRHS(u1, rhs1);

    std::cout << "\nAfter 1 step, element boundary values:\n";
    for (int k = 0; k < mesh.numElements(); ++k)
        std::cout << "  k=" << k
                  << "  u_left="  << u1.elementPtr(k)[0]
                  << "  u_right=" << u1.elementPtr(k)[Np-1] << "\n";

    // Check if discontinuities have appeared
    std::cout << "\nJumps at internal faces after 1 step:\n";
    for (int f = 1; f < mesh.numElements(); ++f) {
        Real jump = u1.elementPtr(f)[0] - u1.elementPtr(f-1)[Np-1];
        std::cout << "  face " << f << ": jump = " << jump << "\n";
    }
//#endif

    auto callback = [&](Real t, const DG::FieldVector<Real>& uh, int step) -> bool {
        if (step % 20 == 0) {
            DG::writeTimeSeriesVTK(mesh, ref, uh, "output", frame++, t);
            std::cout << "t = " << t << "  frame " << frame - 1 << "\n";
        }
        return true;
    };

    DG::Integrator<Real> rk(op, mesh.numElements(), Np, callback);
    rk.integrate(u, 0.0, Tf, dt);

    // Write final state
    DG::writeTimeSeriesVTK(mesh, ref, u, "output", frame++, Tf);
    std::cout << "Done. Written " << frame << " frames.\n";
}
