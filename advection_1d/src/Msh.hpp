#pragma once

#include <TNL/Containers/Array.h>
#include <TNL/Containers/ArrayView.h>
#include <cmath>

// Define TNL namespaces for brevity
using namespace TNL;
using namespace TNL::Containers;

template< typename Real = double, typename Device = Devices::Host, typename Index = int >
struct Mesh2D {
    // Basic Grid Dimensions
    Index K;          // Number of elements
    Index N;          // Polynomial order
    Index Np;         // Number of volume nodes per element
    Index Nfp;        // Number of face nodes per face
    Index Nfaces = 3; // Number of faces per element (triangles)
    
    // Constant tolerance for node matching
    const Real NODETOL = 1e-12; 

    // Physical node coordinates (Size: Np * K)
    Array<Real, Device, Index> x, y;

    // Element-to-Vertex connectivity (Size: K * 3)
    Array<Index, Device, Index> EToV;

    // Metric Elements for local mappings (Size: Np * K)
    Array<Real, Device, Index> rx, sx, ry, sy, J;

    // Surface Normals and Surface Jacobians (Size: Nfaces * Nfp * K)
    Array<Real, Device, Index> nx, ny, sJ, Fscale;

    // Connectivity Arrays (Size: K * Nfaces)
    Array<Index, Device, Index> EToE, EToF;

    // Trace Node Maps (Interior and Exterior traces)
    Array<Index, Device, Index> vmapM, vmapP, mapB, vmapB;

    /**
     * @brief Computes the metric elements for the local mappings of the elements.
     * Corresponds to GeometricFactors2D.m
     */
    void GeometricFactors2D(const Array<Real, Device, Index>& Dr, 
                            const Array<Real, Device, Index>& Ds) 
    {
        Index totalNodes = Np * K;
        rx.setSize(totalNodes); sx.setSize(totalNodes);
        ry.setSize(totalNodes); sy.setSize(totalNodes);
        J.setSize(totalNodes);

        // Assuming Dr and Ds are applied as matrix multiplications over x and y
        // xr = Dr*x; xs = Ds*x; yr = Dr*y; ys = Ds*y; 
        // J = -xs.*yr + xr.*ys;
        // rx =  ys./J; sx = -yr./J; 
        // ry = -xs./J; sy =  xr./J;
        
        // (Implementation of the matrix-vector multiplication using TNL goes here)
    }

    /**
     * @brief Computes outward pointing normals at element faces and surface Jacobians.
     * Corresponds to Normals2D.m
     */
    void Normals2D(const Array<Index, Device, Index>& Fmask) 
    {
        Index totalFaceNodes = Nfaces * Nfp * K;
        nx.setSize(totalFaceNodes); 
        ny.setSize(totalFaceNodes);
        sJ.setSize(totalFaceNodes);
        Fscale.setSize(totalFaceNodes);

        // Logic interpolates volume geometric factors (xr, xs, yr, ys) to face nodes 
        // using Fmask, computes normals based on face definitions:
        // Face 1: nx = fyr, ny = -fxr
        // Face 2: nx = fys - fyr, ny = -fxs + fxr
        // Face 3: nx = -fys, ny = fxs
        // 
        // sJ = sqrt(nx.*nx + ny.*ny);
        // nx = nx./sJ; ny = ny./sJ;
        // Fscale = sJ ./ J(Fmask)
    }

    /**
     * @brief Builds global connectivity arrays for grid based on EToV array.
     * Corresponds to Connect2D.m
     */
    void Connect2D() 
    {
        EToE.setSize(K * Nfaces);
        EToF.setSize(K * Nfaces);

        // In the sources, this is built by creating a sparse face-to-vertex 
        // matrix (SpFToV) and computing SpFToF = SpFToV * SpFToV' - 2*I. 
        // Matches are found where SpFToF == 2.
        // For C++, this translates to hashing the vertex pairs of each face 
        // and matching the faces that share the exact same vertices.
    }

    /**
     * @brief Builds connectivity and boundary tables (vmapM, vmapP).
     * Corresponds to BuildMaps2D.m
     */
    void BuildMaps2D(const Array<Index, Device, Index>& Fmask) 
    {
        // Allocates and builds vmapM (minus/interior traces) by extracting
        // the global node IDs corresponding to the boundary faces using Fmask.
        // 
        // Allocates and builds vmapP (plus/exterior traces) by mapping 
        // neighbor face nodes through the EToE and EToF arrays, and verifying
        // distance with distance matrix D = (x1 - x2)^2 + (y1 - y2)^2 < NODETOL
        
        // Identifies mapB (boundary nodes) where vmapM == vmapP
    }
};
