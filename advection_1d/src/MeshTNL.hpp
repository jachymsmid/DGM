#pragma once

#include <TNL/Meshes/Mesh.h>
#include <TNL/Meshes/MeshBuilder.h>
#include <TNL/Meshes/DefaultConfig.h>
#include <TNL/Meshes/Topologies/Edge.h>
#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>
#include <TNL/Devices/Host.h>
#include <TNL/Devices/Cuda.h>
#include <TNL/Algorithms/parallelFor.h>
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


#include <iostream>


// namespace MeshUtils {
// 
//   using CellTopology = TNL::Meshes::Topologies::Edge;
// 
//   struct MyMeshConfig : public TNL::Meshes::DefaultConfig< CellTopology > {
//     static constexpr bool subentityStorage( int entityDimension, int subentityDimension ) {
//       return subentityDimension == 0 && entityDimension == CellTopology::dimension;
//     }
//   };
// 
//   using MeshType = TNL::Meshes::Mesh< MyMeshConfig >;
// 
//   inline MeshType create1DMesh(int numPoints = 20, double spacing = 0.1) {
//       
//     TNL::Meshes::MeshBuilder< MeshType > builder;
//     const int numCells = numPoints - 1;
//     
//     // FIX 1: Use the new unified allocation method
//     builder.setEntitiesCount(numPoints, numCells);
//     
//     for(int i = 0; i < numPoints; ++i) {
//       builder.setPoint(i, { static_cast<double>(i) * spacing });
//     }
//     
//     for(int i = 0; i < numCells; ++i) {
//       // FIX 2: Remove the '&' because getCellSeed now returns by value
//       auto seed = builder.getCellSeed(i);
//       seed.setCornerId(0, i);
//       seed.setCornerId(1, i + 1);
//     }
//     
//     MeshType mesh;
//     
//     // FIX 3: Call build() directly since it now returns void
//     builder.build(mesh);
// 
//     return mesh;
//   }
// } // namespace MeshUtils


// reference element
template< typename Real, typename Device = TNL::Devices::Host, typename Index = int >
class ReferenceElement1D {
public:
    Index N, Np;

    TNL::Containers::Array<Real, Device, Index> r;
    TNL::Containers::Array<Real, Device, Index> Dr;
    TNL::Containers::Array<Real, Device, Index> LIFT;
    
    // V is usually only needed during setup on the Host
    TNL::Containers::Array<Real, TNL::Devices::Host, Index> V;
};

template< typename Real, typename Device = TNL::Devices::Host, typename Index = int >
class DGMesh1D {
public:
    Index K, Np;

    // Physical coordinates for all nodes. Size: K * Np
    // Flattened 1D array for GPU performance (SoA layout)
    TNL::Containers::Array<Real, Device, Index> x;
    
    // ... connectivity and geometric factors below ...
    // Size: K * 2 (since Nfaces = 2 in 1D)
    TNL::Containers::Array<Index, Device, Index> EToE;
    TNL::Containers::Array<Index, Device, Index> EToF;

    // Holds the global node indices for the boundary points of all elements.
    // Size: K * 2
    TNL::Containers::Array<Index, Device, Index> vmapM;
    TNL::Containers::Array<Index, Device, Index> vmapP;

    // Size: K * Np (or just K if you optimize for affine-only meshes)
    TNL::Containers::Array<Real, Device, Index> rx;
    TNL::Containers::Array<Real, Device, Index> J;

    // Size: K * 2
    TNL::Containers::Array<Real, Device, Index> nx;

    void buildMesh(const std::vector<Real> &VX, Index polynomialOrder)
    {

    }
};

template< typename Real, typename Index = int >
void buildPeriodicConnectivity1D(
    Index K,
    Index N,
    TNL::Containers::Array<Index, TNL::Devices::Host, Index>& EToE,
    TNL::Containers::Array<Index, TNL::Devices::Host, Index>& EToF,
    TNL::Containers::Array<Index, TNL::Devices::Host, Index>& vmapM,
    TNL::Containers::Array<Index, TNL::Devices::Host, Index>& vmapP)
{
    Index Np = N + 1;
    Index Nfaces = 2; // Left face (0) and Right face (1)

    // 1. Resize arrays for K elements
    EToE.setSize(K * Nfaces);
    EToF.setSize(K * Nfaces);
    vmapM.setSize(K * Nfaces);
    vmapP.setSize(K * Nfaces);

    // 2. Build Element-to-Element (EToE) and Element-to-Face (EToF) maps
    for (Index k = 0; k < K; ++k) {
        Index leftFaceIdx = k * Nfaces + 0;
        Index rightFaceIdx = k * Nfaces + 1;

        // --- Left Face (Face 0) ---
        // Connects to the Right Face (Face 1) of the left neighbor (k-1)
        // Periodic boundary: if k==0, left neighbor is K-1
        EToE[leftFaceIdx] = (k == 0) ? (K - 1) : (k - 1);
        EToF[leftFaceIdx] = 1;

        // --- Right Face (Face 1) ---
        // Connects to the Left Face (Face 0) of the right neighbor (k+1)
        // Periodic boundary: if k==K-1, right neighbor is 0
        EToE[rightFaceIdx] = (k == K - 1) ? 0 : (k + 1);
        EToF[rightFaceIdx] = 0;
    }

    // 3. Build Node Maps (vmapM and vmapP)
    // These hold the 1D global array indices of the boundary nodes.
    for (Index k = 0; k < K; ++k) {
        Index leftFaceIdx = k * Nfaces + 0;
        Index rightFaceIdx = k * Nfaces + 1;

        // --- vmapM: Interior trace nodes ---
        // The global index of the first node (left) and last node (right) in element k
        vmapM[leftFaceIdx]  = k * Np + 0;
        vmapM[rightFaceIdx] = k * Np + N;

        // --- vmapP: Exterior trace nodes (from neighbors) ---
        // Look up the neighbor element and face we connected to above
        Index leftNeighbor  = EToE[leftFaceIdx];
        Index leftNeighFace = EToF[leftFaceIdx]; // This will be 1 (right face)

        Index rightNeighbor  = EToE[rightFaceIdx];
        Index rightNeighFace = EToF[rightFaceIdx]; // This will be 0 (left face)

        // The node index on the neighbor's right face is their last node (Np-1, or N)
        vmapP[leftFaceIdx]  = leftNeighbor * Np + N;

        // The node index on the neighbor's left face is their first node (0)
        vmapP[rightFaceIdx] = rightNeighbor * Np + 0;
    }

    std::cout << "Successfully built 1D periodic connectivity for " << K << " elements." << std::endl;
}

template< typename Real, typename Index = int >
void buildGeometricFactors1D(
    Index K, 
    Index N, 
    const TNL::Containers::Array<Real, TNL::Devices::Host, Index>& r,
    const std::vector<Real>& VX, // Array of K+1 element boundary vertices
    TNL::Containers::Array<Real, TNL::Devices::Host, Index>& x,
    TNL::Containers::Array<Real, TNL::Devices::Host, Index>& rx,
    TNL::Containers::Array<Real, TNL::Devices::Host, Index>& J,
    TNL::Containers::Array<Real, TNL::Devices::Host, Index>& nx) 
{
    Index Np = N + 1;
    Index Nfaces = 2; // 1D elements have 2 faces

    // 1. Allocate arrays
    // x, rx, J are evaluated at every volume node
    x.setSize(K * Np);
    rx.setSize(K * Np);
    J.setSize(K * Np);
    
    // nx is evaluated only at the trace nodes (faces)
    nx.setSize(K * Nfaces);

    // 2. Loop over each element to calculate coordinates and metrics
    for (Index k = 0; k < K; ++k) {
        // Get the left and right vertices of element k
        Real x_L = VX[k];
        Real x_R = VX[k + 1];

        // The Jacobian is constant for an affine 1D element
        Real element_J = (x_R - x_L) / 2.0;
        Real element_rx = 1.0 / element_J;

        // Loop over the reference nodes to generate physical nodes
        for (Index i = 0; i < Np; ++i) {
            Index idx = k * Np + i; // 1D flattened index

            // Affine mapping from reference 'r' to physical 'x'
            x[idx] = x_L + 0.5 * (1.0 + r[i]) * (x_R - x_L);
            
            // Populate the metrics (duplicating the constant value at each node)
            J[idx] = element_J;
            rx[idx] = element_rx;
        }

        // 3. Set the surface normals
        Index leftFaceIdx = k * Nfaces + 0;
        Index rightFaceIdx = k * Nfaces + 1;

        // Outward facing normal for the left face is always -1
        nx[leftFaceIdx] = -1.0; 
        
        // Outward facing normal for the right face is always +1
        nx[rightFaceIdx] = 1.0; 
    }

    std::cout << "Successfully built geometric factors for " << K << " elements." << std::endl;
}

// Assuming you have your Host mesh and Reference element already built:
// DGMesh1D<double, TNL::Devices::Host> hostMesh;
// ReferenceElement1D<double, TNL::Devices::Host> hostRef;

template< typename Real, typename Index = int >
void copyMeshToGPU(
    const DGMesh1D<Real, TNL::Devices::Host>& hostMesh,
    const ReferenceElement1D<Real, TNL::Devices::Host>& hostRef,
    DGMesh1D<Real, TNL::Devices::Cuda>& deviceMesh,
    ReferenceElement1D<Real, TNL::Devices::Cuda>& deviceRef)
{
    // 1. Copy Reference Element matrices
    deviceRef.N = hostRef.N;
    deviceRef.Np = hostRef.Np;
    deviceRef.Dr = hostRef.Dr; // Deep copy Host -> Device
    
    // 2. Copy Mesh Dimensions and Factors
    deviceMesh.K = hostMesh.K;
    deviceMesh.Np = hostMesh.Np;
    deviceMesh.rx = hostMesh.rx; // Deep copy Host -> Device
    
    // ... copy x, J, EToE, vmapM, etc ...
}

template< typename Real, typename Device = TNL::Devices::Cuda, typename Index = int >
void computeDerivative1D(
    Index K,
    Index Np,
    const TNL::Containers::Array<Real, Device, Index>& u,   // Input field
    const TNL::Containers::Array<Real, Device, Index>& rx,  // Metric
    const TNL::Containers::Array<Real, Device, Index>& Dr,  // Ref Derivative Matrix
    TNL::Containers::Array<Real, Device, Index>& ux)        // Output derivative field
{
    // 1. Get lightweight views for the GPU lambda
    auto view_u  = u.getConstView();
    auto view_rx = rx.getConstView();
    auto view_Dr = Dr.getConstView();
    auto view_ux = ux.getView();

    // 2. Launch the kernel over all nodes (K * Np threads)
    TNL::Algorithms::parallelFor<Device>(0, K * Np,
        [=] __device__ (Index idx) {
            // Figure out which element (k) and which local node (i) this thread is handling
            Index k = idx / Np;
            Index i = idx % Np;

            Real du_dr = 0.0;

            // Perform the matrix-vector multiplication: Dr * u
            // Dr is size (Np x Np), u for this element is size (Np)
            for(Index j = 0; j < Np; ++j) {
                // Flattened index for the Dr matrix: row i, column j
                Index dr_idx = i * Np + j;
                
                // Flattened index for the u vector: element k, node j
                Index u_idx = k * Np + j;
                
                du_dr += view_Dr[dr_idx] * view_u[u_idx];
            }

            // Apply the chain rule metric and store the result
            view_ux[idx] = view_rx[idx] * du_dr;
        }
    );
}

template< typename Real, typename Device = TNL::Devices::Cuda, typename Index = int >
void computeAdvectionRHS1D(
    Index K, Index Np, Real a_speed,
    const TNL::Containers::Array<Real, Device, Index>& u,
    const TNL::Containers::Array<Real, Device, Index>& rx,
    const TNL::Containers::Array<Real, Device, Index>& nx,
    const TNL::Containers::Array<Index, Device, Index>& vmapM,
    const TNL::Containers::Array<Index, Device, Index>& vmapP,
    const TNL::Containers::Array<Real, Device, Index>& Dr,
    const TNL::Containers::Array<Real, Device, Index>& LIFT,
    TNL::Containers::Array<Real, Device, Index>& rhs) // Output: du/dt
{
    Index Nfaces = 2;
    
    // We need a temporary array to store the flux jumps at the faces.
    // Size: K * Nfaces
    TNL::Containers::Array<Real, Device, Index> fluxJumps(K * Nfaces);
    
    // --- LIGHTWEIGHT VIEWS FOR GPU LAMBDAS ---
    auto view_u = u.getConstView();
    auto view_rx = rx.getConstView();
    auto view_nx = nx.getConstView();
    auto view_vmapM = vmapM.getConstView();
    auto view_vmapP = vmapP.getConstView();
    auto view_Dr = Dr.getConstView();
    auto view_LIFT = LIFT.getConstView();
    auto view_fluxJumps = fluxJumps.getView();
    auto view_rhs = rhs.getView();

    // =====================================================================
    // KERNEL 1: Evaluate Surface Fluxes (Launch over K * Nfaces)
    // =====================================================================
    TNL::Algorithms::parallelFor<Device>(0, K * Nfaces,
        [=] __device__ (Index face_idx) {
            // Get global 1D indices for the minus and plus sides
            Index idM = view_vmapM[face_idx];
            Index idP = view_vmapP[face_idx];

            // Extract the field values
            Real uM = view_u[idM];
            Real uP = view_u[idP];
            Real n_x = view_nx[face_idx];

            // Calculate the difference: delta u = u^+ - u^-
            Real du = uP - uM;

            // Lax-Friedrichs / Upwind Flux Jump formula
            // jump = (a/2) * n_x * du - (|a|/2) * du
            Real jump = 0.5 * a_speed * n_x * du - 0.5 * std::abs(a_speed) * du;

            // Store the computed jump
            view_fluxJumps[face_idx] = jump;
        }
    );

    // =====================================================================
    // KERNEL 2: Volume Derivative & LIFT Operator (Launch over K * Np)
    // =====================================================================
    TNL::Algorithms::parallelFor<Device>(0, K * Np,
        [=] __device__ (Index idx) {
            Index k = idx / Np; // Element index
            Index i = idx % Np; // Local node index inside element

            // 1. Compute the volume derivative contribution (Dr * u)
            Real du_dr = 0.0;
            for(Index j = 0; j < Np; ++j) {
                Index dr_idx = i * Np + j;
                Index u_idx = k * Np + j;
                du_dr += view_Dr[dr_idx] * view_u[u_idx];
            }
            // Chain rule for physical volume derivative
            Real d_dx = view_rx[idx] * du_dr;

            // 2. Compute the LIFT contribution from the faces
            Real lift_contrib = 0.0;
            for(Index f = 0; f < Nfaces; ++f) {
                // LIFT matrix size is (Np x Nfaces)
                Index lift_idx = i * Nfaces + f;
                
                // Fetch the flux jump for this element's face
                Index face_idx = k * Nfaces + f;
                
                lift_contrib += view_LIFT[lift_idx] * view_fluxJumps[face_idx];
            }
            
            // Apply the surface scale factor (which is also rx in 1D)
            lift_contrib = view_rx[idx] * lift_contrib;

            // 3. Assemble the final RHS for the Advection Equation
            // du/dt = -a * du/dx + LIFT(flux_jump)
            view_rhs[idx] = -a_speed * d_dx + lift_contrib;
        }
    );
}

template< typename Real, typename Device = TNL::Devices::Cuda, typename Index = int >
void solveAdvection1D(
    Index K, Index Np, 
    Real finalTime, Real dt, Real a_speed,
    TNL::Containers::Array<Real, Device, Index>& u, // In/Out: The solution field
    const TNL::Containers::Array<Real, Device, Index>& rx,
    const TNL::Containers::Array<Real, Device, Index>& nx,
    const TNL::Containers::Array<Index, Device, Index>& vmapM,
    const TNL::Containers::Array<Index, Device, Index>& vmapP,
    const TNL::Containers::Array<Real, Device, Index>& Dr,
    const TNL::Containers::Array<Real, Device, Index>& LIFT)
{
    Index totalNodes = K * Np;

    // 1. Allocate Low-Storage RK arrays
    TNL::Containers::Array<Real, Device, Index> rhs(totalNodes);
    TNL::Containers::Array<Real, Device, Index> res(totalNodes);
    res.setValue(0.0); // Initialize residual to zero

    // 2. Carpenter & Kennedy 5-stage, 4th-order LSRK coefficients
    const Real rka[5] = {0.0,
                        -567301805773.0 / 1357537059087.0,
                        -2404267990393.0 / 2016746695238.0,
                        -3550918686646.0 / 2091501179385.0,
                        -3270041069962.0 / 2362476167839.0};
                        
    const Real rkb[5] = { 1432997174477.0 / 9575080441755.0,
                          5161836677717.0 / 13612068292357.0,
                          1720146321549.0 / 2090206949498.0,
                          3134564353537.0 / 4481467310338.0,
                          2277821191437.0 / 14882151754819.0};

    // Lightweight views for the update kernel
    auto view_u   = u.getView();
    auto view_res = res.getView();
    auto view_rhs = rhs.getConstView();

    Real time = 0.0;
    int step = 0;

    std::cout << "Starting LSRK4 time integration..." << std::endl;

    // =====================================================================
    // MAIN TIME LOOP
    // =====================================================================
    while (time < finalTime) {
        // Adjust dt for the very last step to hit finalTime exactly
        if (time + dt > finalTime) {
            dt = finalTime - time;
        }

        // Loop over the 5 RK stages
        for (int s = 0; s < 5; ++s) {
            
            // 1. Evaluate the Right-Hand Side (calls the kernel we wrote previously)
            computeAdvectionRHS1D<Real, Device, Index>(
                K, Np, a_speed, u, rx, nx, vmapM, vmapP, Dr, LIFT, rhs
            );

            // Fetch the coefficients for this stage
            Real rka_s = rka[s];
            Real rkb_s = rkb[s];

            // 2. Fused GPU Update Kernel
            // Launch one thread per node to update the residual and solution fields
            TNL::Algorithms::parallelFor<Device>(0, totalNodes,
                [=] __device__ (Index idx) {
                    // res = rka * res + dt * rhs
                    view_res[idx] = rka_s * view_res[idx] + dt * view_rhs[idx];

                    // u = u + rkb * res
                    view_u[idx] = view_u[idx] + rkb_s * view_res[idx];
                }
            );
        }

        time += dt;
        step++;

        // Print progress every 100 steps
        if (step % 100 == 0) {
            std::cout << "Step: " << step << " | Time: " << time << std::endl;
        }
    }

    std::cout << "Integration complete. Final time: " << time << std::endl;
}
