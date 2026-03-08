#include <iostream>
#include "src/Mesh.hpp"

int main() {

    // 1. Get the mesh from our utility header
    // We'll create a small mesh: 5 points, spacing of 0.5
    auto mesh = MeshUtils::create1DMesh(5, 0.5);

    // 2. Query basic mesh properties
    // <0> refers to vertices (0D), <1> refers to cells/edges (1D)
    const int numVertices = mesh.getEntitiesCount<0>();
    const int numCells = mesh.getEntitiesCount<1>();

    std::cout << "--- Mesh Info ---" << std::endl;
    std::cout << "Vertices: " << numVertices << std::endl;
    std::cout << "Cells: " << numCells << "\n" << std::endl;

    // 3. Iterate over vertices to read their coordinates
    std::cout << "--- Vertex Coordinates ---" << std::endl;
    for (int v = 0; v < numVertices; ++v) {
        auto point = mesh.getPoint(v);
        std::cout << "Vertex " << v << " is at x = " << point.x() << std::endl;
    }
    std::cout << "\n";

    // 4. Iterate over cells (edges) to work with topology
    std::cout << "--- Cell Processing ---" << std::endl;
    for (int c = 0; c < numCells; ++c) {
        // Grab the specific cell entity
        auto cell = mesh.getEntity<1>(c);
        
        // A 1D cell (edge) has 2 corners (subentities of dimension 0)
        // Let's get the global IDs of those two vertices
        int leftVertexID = cell.template getSubentityIndex<0>(0);
        int rightVertexID = cell.template getSubentityIndex<0>(1);

        // Fetch their exact spatial coordinates
        double xLeft = mesh.getPoint(leftVertexID).x();
        double xRight = mesh.getPoint(rightVertexID).x();

        // Perform a calculation: find the cell center
        double center = (xLeft + xRight) / 2.0;

        std::cout << "Cell " << c << " connects Vertex " << leftVertexID 
                  << " and Vertex " << rightVertexID 
                  << " | Center is at x = " << center << std::endl;
    }

    return 0;
}
