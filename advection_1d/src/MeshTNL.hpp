#pragma once

#include <iostream>
#include <TNL/Meshes/Mesh.h>
#include <TNL/Meshes/MeshBuilder.h>
#include <TNL/Meshes/DefaultConfig.h>
#include <TNL/Meshes/Topologies/Edge.h>

namespace MeshUtils {

    using CellTopology = TNL::Meshes::Topologies::Edge;

    struct MyMeshConfig : public TNL::Meshes::DefaultConfig< CellTopology > {
        static constexpr bool subentityStorage( int entityDimension, int subentityDimension ) {
            return subentityDimension == 0 && entityDimension == CellTopology::dimension;
        }
    };

    using MeshType = TNL::Meshes::Mesh< MyMeshConfig >;

    inline MeshType create1DMesh(int numPoints = 20, double spacing = 0.1) {
        
        TNL::Meshes::MeshBuilder< MeshType > builder;
        const int numCells = numPoints - 1; 
        
        // FIX 1: Use the new unified allocation method
        builder.setEntitiesCount(numPoints, numCells);
        
        for(int i = 0; i < numPoints; ++i) {
            builder.setPoint(i, { static_cast<double>(i) * spacing }); 
        }
        
        for(int i = 0; i < numCells; ++i) {
            // FIX 2: Remove the '&' because getCellSeed now returns by value
            auto seed = builder.getCellSeed(i);
            seed.setCornerId(0, i);
            seed.setCornerId(1, i + 1);
        }
        
        MeshType mesh;
        
        // FIX 3: Call build() directly since it now returns void
        builder.build(mesh);

        return mesh;
    }
} // namespace MeshUtils
