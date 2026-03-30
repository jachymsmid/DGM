#pragma once

#include <TNL/Meshes/DefaultConfig.h>
#include <TNL/Meshes/Topologies/Edge.h>
#include <TNL/Meshes/Mesh.h>
#include <TNL/Meshes/MeshBuilder.h>
#include <TNL/Containers/StaticArray.h>

namespace DG {

// inherits from DefaultConfig
template
<
  class Real = double,
  class GlobalIndex = int,
  typename LocalIndex = short
>
struct MeshConfig
    : TNL::Meshes::DefaultConfig<TNL::Meshes::Topologies::Edge, 1, Real, GlobalIndex, LocalIndex>
{
  // dont really understand this ??
  static constexpr bool subentityStorageEnabled(int entityDim, int subentityDim)
  {
      return true;
  }
};

// alias for my mesh type
template
<
  class Real = double,
  class Device = TNL::Devices::Host,
  class Index = int
>
using TNLMesh = TNL::Meshes::Mesh<MeshConfig<Real, Index, short>, Device>;

} // namespace DG
