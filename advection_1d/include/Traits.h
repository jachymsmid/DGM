#pragma once

#include <TNL/Devices/Host.h>
#include <TNL/Math.h>
#include <TNL/Algorithms/parallelFor.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Meshes/Mesh.h>
#include <TNL/Meshes/Readers/VTKReader.h>
#include <TNL/Meshes/Readers/VTUReader.h>
#include <TNL/Meshes/TypeResolver/resolveMeshType.h>
#include <TNL/Meshes/Writers/VTKWriter.h>
#include <TNL/Meshes/Writers/VTUWriter.h>

using namespace TNL::Meshes;
using namespace TNL::Containers;
using namespace TNL::Algorithms;

// ------------------------------------------------------------------------------------------------------------------ //

struct Traits
{
   /**
    * which device to launch the code on
    */
   using Device = TNL::Devices::Host;

   /**
    * types used for configuring the input mesh
    */
   using MeshConfig = DefaultConfig< Topologies::Quadrangle >;
   using MeshType = TNL::Meshes::Mesh< MeshConfig, Device >;
   using HostMeshType = TNL::Meshes::Mesh< MeshConfig, TNL::Devices::Host >;
   using ReaderType = Readers::VTKReader;
   using WriterType = Writers::VTKWriter< HostMeshType >;

   /**
    * types for basic data types
    */
   using RealType = float;
   using IndexType = int;
   using LocalIndexType = MeshType::LocalIndexType;
   using GlobalIndexType = MeshType::GlobalIndexType;

};

// ------------------------------------------------------------------------------------------------------------------ //
