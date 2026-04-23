#include "Simulation.h"
#include "Traits.h"

#include <TNL/Config/parseCommandLine.h>

using RealType  = typename Traits::RealType;
using IndexType = typename Traits::IndexType;

// ------------------------------------------------------------------------------------------------------------------ //

int
main( int argc, char* argv[] )
{
   TNL::Config::ConfigDescription config;
   config.addEntry< std::string >(
      "meshPath", "provide a path to the mesh.vtu file", "meshes/mesh_regular_12.vtk" );
   config.addEntry< IndexType >( "N", "Order of polynomial approximation", 4);
   config.addEntry< RealType >( "CFL", "Courant number", 0.5 );
   config.addEntry< RealType >("T", "Final time", 2);

   TNL::Config::ParameterContainer parameters;
   if( ! TNL::Config::parseCommandLine( argc, argv, config, parameters ) ) {
      return EXIT_FAILURE;
   }

   const auto      meshPath = parameters.getParameter< std::string >( "meshPath" );
   const RealType  CFL      = parameters.getParameter< RealType >( "CFL" );
   const IndexType N        = parameters.getParameter<IndexType>("N");

   Simulation< Traits > simulation( meshPath );

   simulation.setInitialCondition();
   simulation.setPeriodicBoundaryCondition();

   std::cout << "Initial conditions set." << '\n';

   simulation.run();
   simulation.exportResults( "results.vtu" );

   return EXIT_SUCCESS;
}
