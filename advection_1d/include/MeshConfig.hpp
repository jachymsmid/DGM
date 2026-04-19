/**
 * @file MeshConfig.hpp
 * @brief TNL mesh configuration and convenient alias for the 1D mesh type.
 *
 * Defines a MeshConfig that inherits from TNL::Meshes::DefaultConfig for
 * 1D Edge topologies and provides the TNLMesh alias used across the codebase.
 */

#pragma once

#include <TNL/Meshes/DefaultConfig.h>
#include <TNL/Meshes/Topologies/Edge.h>
#include <TNL/Meshes/Mesh.h>
#include <TNL/Meshes/MeshBuilder.h>
#include <TNL/Containers/StaticArray.h>

namespace DG {

// inherits from DefaultConfig
/**
 * @struct MeshConfig
 * @brief TNL mesh configuration for 1D Edge topology.
 *
 * Inherits from TNL::Meshes::DefaultConfig and enables subentity
 * storage for the 1D mesh representation used by the DG solver.
 */
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
/**
 * @brief Alias for the TNL mesh type used in the DG codebase.
 */
template
<
  class Real = double,
  class Device = TNL::Devices::Host,
  class Index = int
>
using TNLMesh = TNL::Meshes::Mesh<MeshConfig<Real, Index, short>, Device>;

} // namespace DG
