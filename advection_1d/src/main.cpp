#include "FieldVector.hpp"
#include "IO.hpp"
#include "Mesh.hpp"
#include "NumericalFlux.hpp"
#include "Integrator.hpp"
#include "Operator.hpp"
#include <TNL/Containers/StaticArray.h>
#include <TNL/Containers/Vector.h>
#include <iostream>

using Device  = TNL::Devices::Host;
using Real = double;
using Index = int;

int main(int argc, char* argv[])
{
    const int  K = 12; // number of elements
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

    // -------------------- initial conditions --------------------------------

    // create array views
    auto u_view = u.data().getView();

    auto sin_init = [=] __cuda_callable__ ( const TNL::Containers::StaticArray< 2, int >& i  ) mutable
    {
      Real xL = mesh.leftVertex(i.x());
      Real h = mesh.elementSize(i.x());
      Real r = ref.nodes()[i.y()];
      u_view[ i.y() + i.x() * ref.numDOF() ] = TNL::sin(xL + (r + 1.0) * 0.5 * h);
    };

    auto saw_init = [=] __cuda_callable__ ( const TNL::Containers::StaticArray< 2, int >& idx) mutable
    {
      if ( idx.x() < int(mesh.numElements()/3) || idx.x() > int(2 * mesh.numElements()/3) )
      {
        u_view[ idx.y() + idx.x() * ref.numDOF() ] = 0.0;
      }
      else
      {
        u_view[ idx.y() + idx.x() * ref.numDOF() ] = 1.0;
      }
    };

    auto cone_init = [=] __cuda_callable__ ( const TNL::Containers::StaticArray< 2, int >& idx ) mutable
    {
      Real xL = mesh.leftVertex(idx.x());
      Real h = mesh.elementSize(idx.x());
      Real r = ref.nodes()[idx.y()];

      if ( idx.x() < int(mesh.numElements()/4) )
      {
        u_view[ idx.y() + idx.x() * ref.numDOF() ] = (xL + (r + 1.0) * 0.5 * h)/( h * int(mesh.numElements()/4) );
      }
      else if ( idx.x() >= int(mesh.numElements()/4) && idx.x() < int(mesh.numElements()/2) )
      {
        u_view[ idx.y() + idx.x() * ref.numDOF() ] = - (xL + (r + 1.0) * 0.5 * h)/( h * int(mesh.numElements()/4) ) + 2;
      }
      else
      {
        u_view[ idx.y() + idx.x() * ref.numDOF() ] = 0.0;
      }
    };

    TNL::Containers::StaticArray< 2, int > begin{0, 0};
    // we expect same number of DOF on each elemnt
    TNL::Containers::StaticArray< 2, int > end{mesh.numElements(), ref.numDOF()};

    TNL::Algorithms::parallelFor< Device >(begin, end, cone_init);

    Real r_min = 2.0;
    for (int k = 0; k < mesh.numElements(); k++)
    {
      for (int i = 1; i < ref.numDOF(); i++)
      {
        r_min = std::min(r_min, ref.nodes()[i] - ref.nodes()[i-1]);
      }
    }

    Real x_min = r_min * mesh.minJacobian();

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
    Real dt = DG::ERK<Real>::computeDt(x_min, a, N);

    auto callback = [&](Real t, const DG::FieldVector<Real>& uh, int step) -> bool
    {
      // DG::writeTimeSeriesVTK(mesh, ref, uh, "output/output", frame++, t);
      return true;
    };

    DG::ERK<Real> rk(op.rhsFunction(), mesh.numElements(), ref.numDOF(), callback);
    rk.integrate(u, 0.0, Tf, dt);

    // Write final state
    DG::writeTimeSeriesVTK(mesh, ref, u, "output/output", frame++, Tf);
    std::cout << "Done. Written " << frame << " frames.\n";
}
