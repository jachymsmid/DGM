#include "Solver.hpp"

using Device  = TNL::Devices::Host;
using Real = double;
using Index = int;

int main()
{
   TNL::DGM::Solver<Model> solver;
   solver.init();
   solver.makeSnapshot();

   while ( solver.timeStepping.runSimulation() )
   {
      solver.step();
      solver.applyFilter();
   }

   solver.makeSnapshot();

   // move to config struct (prameter container)
   const int  K = 12; // number of elements
   const int N = 4; // polynomial order of approximation
   const Real a = 1.0; // advection speed
   const Real Tf = 2.0; // final time
   const Real CFL = 0.4;

   /**
    * TODO: put inside Model object
    */
   static auto physical_flux = [&] ( Real u ) -> Real { return a * u; };
   static auto advection_speed = [&] ( Real u ) -> Real { return a; };

   //// -------------------- initial conditions --------------------------------
   //// create array views
   //TNL::DGM::FieldVector u = solver.getData();
   //TNL::DGM::Mesh mesh = solver.getMesh();

   //auto u_view = u.data().getView();

   //auto sin_init = [=] __cuda_callable__ ( const TNL::Containers::StaticArray< 2, int >& i  ) mutable
   //{
   //  Real xL = mesh.leftVertex(i.x());
   //  Real h = mesh.elementSize(i.x());
   //  Real r = ref.nodes()[i.y()];
   //  u_view[ i.y() + i.x() * ref.numDOF() ] = - TNL::sin((xL + (r + 1.0) * 0.5 * h) * PI) + 1.f;
   //};

   //auto shock_init = [=] __cuda_callable__ ( const TNL::Containers::StaticArray< 2, int >& idx ) mutable
   //{
   //  if (idx.x() < int(mesh.numElements()/3) || idx.x() > int(mesh.numElements()/2))
   //  {
   //    u_view[ idx.y() + idx.x() * ref.numDOF() ] = 1.0;
   //  }
   //  else
   //  {
   //    u_view[ idx.y() + idx.x() * ref.numDOF() ] = 0.0;
   //  }
   //};

   //auto rarefaction_init = [=] __cuda_callable__ ( const TNL::Containers::StaticArray< 2, int >& idx ) mutable
   //{
   //  if (idx.x() < int(mesh.numElements()/3) || idx.x() > int(mesh.numElements()/2))
   //  {
   //    u_view[ idx.y() + idx.x() * ref.numDOF() ] = 0.0;
   //  }
   //  else
   //  {
   //    u_view[ idx.y() + idx.x() * ref.numDOF() ] = 1.0;
   //  }
   //};

   //auto saw_init = [=] __cuda_callable__ ( const TNL::Containers::StaticArray< 2, int >& idx) mutable
   //{
   //  if ( idx.x() < int(mesh.numElements()/3) || idx.x() > int(2 * mesh.numElements()/3) )
   //  {
   //    u_view[ idx.y() + idx.x() * ref.numDOF() ] = 0.0;
   //  }
   //  else
   //  {
   //    u_view[ idx.y() + idx.x() * ref.numDOF() ] = 1.0;
   //  }
   //};

   //auto cone_init = [=] __cuda_callable__ ( const TNL::Containers::StaticArray< 2, int >& idx ) mutable
   //{
   //  Real xL = mesh.leftVertex(idx.x());
   //  Real h = mesh.elementSize(idx.x());
   //  Real r = ref.nodes()[idx.y()];

   //  if ( idx.x() < int(mesh.numElements()/4) )
   //  {
   //    u_view[ idx.y() + idx.x() * ref.numDOF() ] = (xL + (r + 1.0) * 0.5 * h)/( h * int(mesh.numElements()/4)) + 2.0;
   //  }
   //  else if ( idx.x() >= int(mesh.numElements()/4) && idx.x() < int(mesh.numElements()/2) )
   //  {
   //    u_view[ idx.y() + idx.x() * ref.numDOF() ] = - (xL + (r + 1.0) * 0.5 * h)/( h * int(mesh.numElements()/4) );
   //  }
   //  else
   //  {
   //    u_view[ idx.y() + idx.x() * ref.numDOF() ] = 0.0;
   //  }
   //};

   //TNL::Containers::StaticArray< 2, int > begin{0, 0};
   //// we expect same number of DOF on each elemnt
   //TNL::Containers::StaticArray< 2, int > end{mesh.numElements(), ref.numDOF()};

   //// 2-dimensional parallel for
   //TNL::Algorithms::parallelFor< Device >(begin, end, saw_init);
}
