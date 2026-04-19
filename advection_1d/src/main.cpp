#include "FieldVector.hpp"
#include "IO.hpp"
#include "Mesh.hpp"
#include "NumericalFlux.hpp"
#include "Integrator.hpp"
#include "Operator.hpp"
#include <TNL/Containers/StaticArray.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Math.h>
#include <iostream>

using Device  = TNL::Devices::Host;
using Real = double;
using Index = int;

int main(int argc, char* argv[])
{
    const int  K = 20; // number of elements
    const int N = 3; // polynomial order of approximation
    const Real a = 1.0; // advection speed
    const Real Tf = 2.0; // final time
    const Real CFL = 0.4;
    const Real PI = TNL::pi;

    // Build mesh: either from a VTK file or uniform
    // DG::Mesh<Real> mesh = (argc > 1)
    //     ? DG::Mesh<Real>::readVTK(argv[1])
    //     : DG::Mesh<Real>::uniform(0.0, 2.0 * M_PI, K);

    // Burger's equation
    auto physical_flux = [&] ( Real u ) -> Real { return 1.0/2.0 * u * u; };
    auto advection_speed = [&] ( Real u ) -> Real { return u; };
    // Linear advection
    // auto physical_flux = [&] ( Real u ) -> Real { return a * u; };
    // auto advection_speed = [&] ( Real u ) -> Real { return a; };

    // construct uniform mesh
    DG::Mesh<Real> mesh = DG::Mesh<Real>::uniform(-1.0, 1.0, K);

    // construct reference element
    DG::ReferenceElement<Real> ref(N);

    // construct numerical flux: UpwindFlux, LaxFriedrichsFlux, GodunovFlux, RoeFlux
    DG::RoeFlux<Real> numerical_flux( advection_speed, physical_flux );

    // construct rhs operator
    DG::Operator<Real> op(mesh, ref, numerical_flux, physical_flux);

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
      u_view[ i.y() + i.x() * ref.numDOF() ] = TNL::sin((xL + (r + 1.0) * 0.5 * h) * PI);
    };

    auto abs_sin_init = [=] __cuda_callable__ ( const TNL::Containers::StaticArray< 2, int >& i  ) mutable
    {
      Real xL = mesh.leftVertex(i.x());
      Real h = mesh.elementSize(i.x());
      Real r = ref.nodes()[i.y()];
      u_view[ i.y() + i.x() * ref.numDOF() ] = TNL::abs(TNL::sin((xL + (r + 1.0) * 0.5 * h) * PI));
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

    // 2-dimensional parallel for
    TNL::Algorithms::parallelFor< Device >(begin, end, abs_sin_init);

    // -------------------------- more setup ----------------------------------
    // find delta x_min for time step computation
    Real r_min = 2.0;
    for (int k = 0; k < mesh.numElements(); k++)
    {
      for (int i = 1; i < ref.numDOF(); i++)
      {
        r_min = std::min(r_min, ref.nodes()[i] - ref.nodes()[i-1]);
      }
    }

    Real x_min = r_min * mesh.minJacobian();

    // auto energy = [&](const DG::FieldVector<Real>& v) {
    //   Real E = 0;
    //   for (int k = 0; k < mesh.numElements(); ++k) {
    //       Real J = mesh.jacobian(k);
    //       const Real* vk = v.elementPtr(k);
    //       for (int i = 0; i < Np; ++i)
    //           E += ref.weights()[i] * J * vk[i] * vk[i];
    //   }
    //   return E;
    // };

    // find max advection speed
    Real max_speed = 0.0;
    for (int i = 0; i < mesh.numElements(); i++)
    {
      for (int j = 0; j < ref.numDOF(); j++)
      {
        max_speed = TNL::argAbsMax(max_speed, advection_speed(u.elementPtr(i)[j]));
      }
    }

    Real dt = DG::ERK<Real>::computeDt(x_min, max_speed, N, CFL);

    // write initial condition
    int frame = 0;
    DG::writeTimeSeriesVTK(mesh, ref, u, "output/output", frame++, Real(0));

    auto callback = [&](Real t, const DG::FieldVector<Real>& uh, int step) -> bool
    {
      DG::writeTimeSeriesVTK(mesh, ref, uh, "output/output", frame++, t);
      return true;
    };

    DG::SSPRK<Real> rk(op.rhsFunction(), mesh.numElements(), ref.numDOF(), callback);
    rk.integrate(u, 0.0, Tf, dt);

    // output
    std::cout << "Starting simulation with: " << std::endl;
    std::cout << "\tK = " << K << std::endl;
    std::cout << "\tN = " << N << std::endl;
    std::cout << "\tx_min = " << x_min << std::endl;
    std::cout << "\tdt = " << dt << std::endl;
    std::cout << "\tmax advection speed = " << max_speed << std::endl;
    DG::writeTimeSeriesVTK(mesh, ref, u, "output/output", frame++, Tf);
    std::cout << "Done. Written " << frame << " frames.\n";
}
