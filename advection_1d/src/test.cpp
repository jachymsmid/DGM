#define NDBUG
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include <TNL/Assert.h>
#include <TNL/Containers/Array.h>
#include <TNL/Devices/Host.h>
#include <TNL/Devices/Cuda.h>
#include <TNL/Meshes/Mesh.h>
#include <TNL/Meshes/MeshBuilder.h>
#include <TNL/Meshes/DefaultConfig.h>
#include <TNL/Meshes/Topologies/Edge.h>
#include <TNL/Containers/ArrayView.h>
#include <TNL/Devices/Host.h>
#include <TNL/Devices/Cuda.h>
#include <TNL/Algorithms/parallelFor.h>

#include "MeshTNL.hpp"

/**
 * @brief Helper to generate a uniform 1D mesh of vertices
 */
std::vector<double> createUniformMesh1D(double xL, double xR, int K) {
    std::vector<double> VX(K + 1);
    double dx = (xR - xL) / K;
    for (int i = 0; i <= K; ++i) {
        VX[i] = xL + i * dx;
    }
    return VX;
}

void writeSolutionToCSV(const std::string& filename, int K, int Np,
                        const TNL::Containers::Array<double, TNL::Devices::Host>& x,
                        const TNL::Containers::Array<double, TNL::Devices::Host>& u)
{
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write the CSV header
    file << "x,u\n";

    // Write the data with high precision
    int totalNodes = K * Np;
    for (int i = 0; i < totalNodes; ++i) {
        file << std::fixed << std::setprecision(8) << x[i] << "," << u[i] << "\n";
    }

    file.close();
    std::cout << "Successfully saved solution to " << filename << std::endl;
}

int main() {
    // --- 1. Simulation Parameters ---
    int K = 20;                 // Number of elements
    int N = 4;                  // Polynomial order
    int Np = N + 1;             // Nodes per element
    double a_speed = 1.0;       // Advection wave speed
    double finalTime = 2.0;     // Run for one full period on domain [-1, 1]

    // Calculate a stable timestep based on the DG CFL condition
    // For DG, dt scales with dx / (speed * N^2)
    double dx = 2.0 / K;
    double dt = 0.25 * dx / (a_speed * N * N);

    std::cout << "Setting up DG Advection 1D..." << std::endl;
    std::cout << "Elements: " << K << " | Order: " << N << " | dt: " << dt << std::endl;

    // --- 2. Host Setup (Reference Element) ---
    TNL::Containers::Array<double, TNL::Devices::Host> host_r(Np);
    TNL::Containers::Array<double, TNL::Devices::Host> host_Dr(Np * Np);
    TNL::Containers::Array<double, TNL::Devices::Host> host_LIFT(Np * 2);

    // *************************************************************************
    // NOTE: In a complete implementation, you would populate host_r (GLL nodes),
    // host_Dr, and host_LIFT here using Vandermonde matrix operations.
    // For testing the pipeline, ensure these arrays are mathematically valid!
    // *************************************************************************

    // --- 3. Host Setup (Physical Mesh) ---
    std::vector<double> VX = createUniformMesh1D(-1.0, 1.0, K);

    TNL::Containers::Array<double, TNL::Devices::Host> host_x(K * Np);
    TNL::Containers::Array<double, TNL::Devices::Host> host_rx(K * Np);
    TNL::Containers::Array<double, TNL::Devices::Host> host_J(K * Np);
    TNL::Containers::Array<double, TNL::Devices::Host> host_nx(K * 2);

    TNL::Containers::Array<int, TNL::Devices::Host> host_EToE(K * 2);
    TNL::Containers::Array<int, TNL::Devices::Host> host_EToF(K * 2);
    TNL::Containers::Array<int, TNL::Devices::Host> host_vmapM(K * 2);
    TNL::Containers::Array<int, TNL::Devices::Host> host_vmapP(K * 2);

    buildPeriodicConnectivity1D<double, int>(K, N, host_EToE, host_EToF, host_vmapM, host_vmapP);
    buildGeometricFactors1D<double, int>(K, N, host_r, VX, host_x, host_rx, host_J, host_nx);

    // --- 4. Initial Condition (Gaussian Pulse) ---
    TNL::Containers::Array<double, TNL::Devices::Host> host_u(K * Np);
    for (int i = 0; i < K * Np; ++i) {
        // Evaluate the Gaussian pulse at each physical coordinate
        host_u[i] = std::exp(-40.0 * host_x[i] * host_x[i]);
    }

    // --- 5. Transfer to Device (GPU) ---
    std::cout << "Transferring data to GPU..." << std::endl;
    // Overloaded assignment operators in TNL automatically trigger cudaMemcpy
    TNL::Containers::Array<double, TNL::Devices::Cuda> dev_x = host_x;
    TNL::Containers::Array<double, TNL::Devices::Cuda> dev_rx = host_rx;
    TNL::Containers::Array<double, TNL::Devices::Cuda> dev_nx = host_nx;
    TNL::Containers::Array<double, TNL::Devices::Cuda> dev_u = host_u;
    TNL::Containers::Array<double, TNL::Devices::Cuda> dev_Dr = host_Dr;
    TNL::Containers::Array<double, TNL::Devices::Cuda> dev_LIFT = host_LIFT;
    
    TNL::Containers::Array<int, TNL::Devices::Cuda> dev_vmapM = host_vmapM;
    TNL::Containers::Array<int, TNL::Devices::Cuda> dev_vmapP = host_vmapP;

    // --- 6. Run the Solver on the GPU ---
    solveAdvection1D<double, TNL::Devices::Cuda, int>(
        K, Np, finalTime, dt, a_speed, 
        dev_u, dev_rx, dev_nx, dev_vmapM, dev_vmapP, dev_Dr, dev_LIFT
    );

    // --- 7. Bring Results Back to CPU ---
    std::cout << "Retrieving solution from GPU..." << std::endl;
    host_u = dev_u; // Copies the final timestep back to the CPU

    // --- 8. Write to CSV ---
    writeSolutionToCSV("advection_results.csv", K, Np, host_x, host_u);

    std::cout << "Simulation finished successfully!" << std::endl;

    return 0;
}
