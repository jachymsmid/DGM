#include "IO.hpp"
#include "Mesh.hpp"
#include "NumericalFlux.hpp"
#include "RK4Integrator.hpp"
#include <TNL/Containers/Vector.h>
#include <iostream>
#include <cmath>

int main(int argc, char* argv[])
{
    using Real = double;
    const int  K = 10; // number of elements
    const int N = 4; // polynomial order of approximation
    const Real a = 1.0; // advection speed
    const Real Tf = 2.0 * M_PI; // final time

    // Build mesh: either from a VTK file or uniform
    // DG::Mesh<Real> mesh = (argc > 1)
    //     ? DG::Mesh<Real>::readVTK(argv[1])
    //     : DG::Mesh<Real>::uniform(0.0, 2.0 * M_PI, K);

    // construct uniform mesh
    DG::Mesh<Real> mesh = DG::Mesh<Real>::uniform(0.0, 2.0 * M_PI, K);

    // construct reference element
    DG::ReferenceElement<Real> ref(N);

    // construct numerical flux
    // DG::UpwindFlux<Real> flux(a);
    DG::UpwindFlux<Real> numerical_flux(a);

    // construct rhs operator
    DG::Operator<Real> op(mesh, ref, numerical_flux, [&](Real u){ return a * u; });

    // construct field vector
    DG::FieldVector<Real> u(mesh.numElements(), ref.numDOF());

    int Np = ref.numDOF();

    // initial condition u(x,0) = sin(x)
    for (int k = 0; k < mesh.numElements(); ++k) {
      Real xL = mesh.leftVertex(k), h = mesh.elementSize(k);
      Real* uk = u.elementPtr(k);
      for (int i = 0; i < ref.numDOF(); ++i) {
        Real r = ref.nodes()[i];
        // maybe wrap the affine mapping in a member function
        uk[i]  = std::sin(xL + (r + 1.0) * 0.5 * h);
      }
    }

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

    // write initial condition
    int frame = 0;
    DG::writeTimeSeriesVTK(mesh, ref, u, "output/output", frame++, Real(0));

    // Time stepping with output every 20 steps
    Real h_min = mesh.minElementSize();
    std::cout << "min(h) = " << h_min << std::endl;
    Real dt = DG::Integrator<Real>::computeDt(h_min, a, N);
    std::cout << "Time step for intergation:  " << dt << std::endl;
    std::cout << "Energy at t = 0 :" << energy(u) << std::endl;

    auto callback = [&](Real t, const DG::FieldVector<Real>& uh, int step) -> bool
    {
      DG::writeTimeSeriesVTK(mesh, ref, uh, "output/output", frame++, t);

      // --- diagnostics ---
      if (step == 20)
      {
        std::cout << "h_min = " << h_min << "  dt = " << dt << "\n\n";

        TNL::Containers::Vector< Real, TNL::Devices::Host, int > r,w;
        // ref.compute_printGLL(r, w, 10);

        // sum of weights should equal 2.0 (length of [-1,1])
        Real wsum = 0;
        for (int i = 0; i < Np; ++i) wsum += ref.weights()[i];
        std::cout << "sum(weights) = " << wsum << "  (should be 2.0)\n\n";

        // D matrix: D*ones should be zero (derivative of constant = 0)
        std::cout << "D * ones (should all be ~0):\n";
        for (int i = 0; i < Np; ++i) {
            Real s = 0;
            for (int j = 0; j < Np; ++j) s += ref.Dr()(i,j);
            std::cout << "  row " << i << ": " << s << "\n";
        }

        // D matrix: D*r should be ones (derivative of identity = 1)
        std::cout << "D * r (should all be ~1.0):\n";
        for (int i = 0; i < Np; ++i) {
            Real s = 0;
            for (int j = 0; j < Np; ++j) s += ref.Dr()(i,j) * ref.nodes()[j];
            std::cout << "  row " << i << ": " << s << "\n";
        }

        // LIFT: columns should sum to Np (surface-to-volume consistency)
        std::cout << "\nLIFT col 0 (left face):  ";
        Real lsum0 = 0;
        for (int i = 0; i < Np; ++i) { std::cout << ref.LIFT()(i,0) << " "; lsum0 += ref.LIFT()(i,0); }
        std::cout << "\n  sum = " << lsum0 << "\n";

        std::cout << "LIFT col 1 (right face): ";
        Real lsum1 = 0;
        for (int i = 0; i < Np; ++i) { std::cout << ref.LIFT()(i,1) << " "; lsum1 += ref.LIFT()(i,1); }
        std::cout << "\n  sum = " << lsum1 << "\n";

        std::cout << "Energy at t = 20 :" << energy(u) << std::endl;

        std::cout << "\nJumps at internal faces after 20 steps:\n";
        for (int f = 1; f < mesh.numElements(); ++f) {
          Real* ul = u.elementPtr(f);
          Real* ur = u.elementPtr(f-1);
            Real jump = ul[0] - ur[Np-1];
            std::cout << "  face " << f << ": jump = " << jump << "\n";
        }

      }
      return true;
    };

    DG::Integrator<Real> rk(op, mesh.numElements(), ref.numDOF(), callback);
    rk.integrate(u, 0.0, Tf, dt);

    // Write final state
    DG::writeTimeSeriesVTK(mesh, ref, u, "output/output", frame++, Tf);
    std::cout << "Done. Written " << frame << " frames.\n";
}
