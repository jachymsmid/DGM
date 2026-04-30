/**
 * @file IO.hpp
 * @brief VTK input/output helpers for writing DG solutions.
 *
 * Contains utilities to write DG data as legacy ASCII VTK files so that
 * Paraview can visualize element-local polynomial solutions and time
 * series (writeToVTK, writeTimeSeriesVTK).
 */

#pragma once

#include "Mesh.hpp"
#include "ReferenceElement.hpp"
#include "FieldVector.hpp"
#include "PadeLegendre.hpp"
#include <TNL/Meshes/Readers/VTKReader.h>
#include <TNL/Meshes/Writers/VTKWriter.h>
#include <fstream>
#include <string>
#include <stdexcept>
#include <iomanip>

namespace DG {

// ─────────────────────────────────────────────────────────────────────────────
//  writeToVTK
//
//  Writes the DG solution to a legacy ASCII VTK file that Paraview can open.
//
//  Strategy: the DG solution lives on sub-element nodes, not on the mesh
//  vertices. We therefore write a *refined* point cloud: for each element k
//  we emit its Np physical node positions as VTK points, and connect each
//  group of Np points into a VTK_POLY_LINE cell. This gives Paraview a
//  piecewise-polynomial picture of u inside every element, with the expected
//  discontinuities visible at element boundaries.
//
//  Output layout:
//    - K * Np  points  (physical DG node positions)
//    - K       cells   (one VTK_POLY_LINE per element, with Np vertices)
//    - 1 point data array named after `fieldName`
//
//  Usage:
//      DG::writeToVTK(mesh, ref, u, "output.vtk", "u");
// ─────────────────────────────────────────────────────────────────────────────
template
<
  class Real = double,
  class Device = TNL::Devices::Host,
  class Index = int
>
/**
 * @brief Write a DG solution to an ASCII legacy VTK file.
 *
 * The function emits K*Np points and K POLY_LINE cells so Paraview can
 * display the piecewise-polynomial DG solution including element
 * discontinuities. The produced file includes a POINT_DATA scalar array
 * with the field values.
 */
void writeToVTK(const Mesh<Real, Device, Index>& mesh,
                const ReferenceElement<Real, Index>& ref,
                const FieldVector<Real, Device, Index>& u,
                const std::string& filename,
                const std::string& fieldName = "u",
                Real t = Real(0))
{
  const Index K = mesh.numElements();
  const Index Np = ref.numDOF();
  const Index totalPts = K * Np;

  std::ofstream f(filename);
  if (!f)
    throw std::runtime_error("writeToVTK: cannot open '" + filename + "'");

  f << std::scientific << std::setprecision(10);

  // ── VTK header ──────────────────────────────────────────────────────────
  f << "# vtk DataFile Version 3.0\n";
  f << "DG1D solution: " << fieldName << " t=" << t << "\n";
  f << "ASCII\n";
  f << "DATASET UNSTRUCTURED_GRID\n";

  // ── Points ──────────────────────────────────────────────────────────────
  // VTK always needs 3D coordinates; pad y=0, z=0 for 1D.
  f << "POINTS " << totalPts << " double\n";
  for (Index k = 0; k < K; ++k)
  {
    const Real xL = mesh.leftVertex(k);
    const Real h  = mesh.elementSize(k);
    for (Index i = 0; i < Np; ++i)
    {
      // map reference node r_i in [-1,1] to physical coordinate
      const Real r = ref.nodes()[i];
      const Real x = xL + (r + Real(1)) * Real(0.5) * h;
      f << x << " 0.0 0.0\n";
    }
  }

  // ── Cells ────────────────────────────────────────────────────────────────
  // Each cell is a VTK_POLY_LINE with Np vertices.
  // Cell list entry: [count, v0, v1, ..., v_{Np-1}]
  const Index cellListSize = K * (Np + 1);  // Np vertices + 1 count per cell
  f << "\nCELLS " << K << " " << cellListSize << "\n";
  for (Index k = 0; k < K; ++k)
  {
    f << Np;
    for (Index i = 0; i < Np; ++i)
    {
      f << " " << k * Np + i;
    }
    f << "\n";
  }

  // VTK_POLY_LINE = cell type 4
  f << "\nCELL_TYPES " << K << "\n";
  for (Index k = 0; k < K; ++k)
  {
    f << "4\n";
  }

  // ── Point data ───────────────────────────────────────────────────────────
  f << "\nPOINT_DATA " << totalPts << "\n";
  f << "SCALARS " << fieldName << " double 1\n";
  f << "LOOKUP_TABLE default\n";
  for (Index k = 0; k < K; ++k)
  {
    const Real* uk = u.elementPtr(k);
    for (Index i = 0; i < Np; ++i)
    {
      f << uk[i] << "\n";
    }
  }
}

// ─────────────────────────────────────────────────────────────────────────────
//  writeTimeSeriesVTK
//
//  Convenience wrapper: writes a numbered sequence of VTK files suitable for
//  loading as a time series in Paraview.
//
//  Call it from your RK4 callback:
//
//    int frame = 0;
//    auto cb = [&](Real t, const Field& u, int step) -> bool {
//        if (step % 10 == 0)
//            DG::writeTimeSeriesVTK(mesh, ref, u, "output", frame++, t);
//        return true;
//    };
//
//  Produces: output_000000.vtk, output_000001.vtk, ...
//  Load them in Paraview with File > Open, select all, apply "Group Datasets".
// ─────────────────────────────────────────────────────────────────────────────
template
<
  class Real,
  class Device,
  class Index
>
/**
 * @brief Convenience wrapper to write numbered VTK frames for a timeseries.
 *
 * Produces files named basename_XXXXXX.vtk suitable for loading as a
 * grouped time series in Paraview.
 */
void writeTimeSeriesVTK(const Mesh<Real, Device, Index>& mesh,
                        const ReferenceElement<Real, Index>& ref,
                        const FieldVector<Real, Device, Index>& u,
                        const std::string& basename,
                        int frameIndex,
                        Real t,
                        const std::string& fieldName = "u")
{
    // Zero-pad frame index to 6 digits: output_000042.vtk
    std::ostringstream oss;
    oss << basename << "_" << std::setfill('0') << std::setw(6) << frameIndex << ".vtk";
    writeToVTK(mesh, ref, u, oss.str(), fieldName, t);
}

// ─────────────────────────────────────────────────────────────────────────────
//  writePadeVTK
//
//  Writes the Padé–Legendre reconstruction to a legacy ASCII VTK file.
//
//  Strategy: for each element k, `refineFactor * Np` equidistant reference
//  points are generated on [-1,1], the Padé approximant is evaluated there,
//  and the points are connected into a VTK_POLY_LINE cell.  This gives a
//  smooth, high-resolution visualization of the rational approximant,
//  independent of the (irregularly-spaced) GLL nodes used by the DG solver.
//
//  Usage:
//      DG::writePadeVTK(mesh, ref, u, pade_solver, "output/pade.vtk", "u_pade", t, 4);
//
//  Parameters:
//    mesh          — the DG mesh
//    ref           — reference element (provides Np)
//    u             — DG solution FieldVector (nodal values at GLL nodes)
//    solver        — a constructed PadeLegendreSolver (same L,M,N as desired)
//    filename      — output file path
//    fieldName     — scalar field name written to the VTK header
//    t             — simulation time (written to the VTK comment line)
//    refineFactor  — number of equidistant sub-points per GLL node (default 4)
// ─────────────────────────────────────────────────────────────────────────────
template
<
  class Real   = double,
  class Device = TNL::Devices::Host,
  class Index  = int
>
/**
 * @brief Write the Padé–Legendre reconstruction at a finer equidistant grid.
 *
 * For each element, `refineFactor * Np` uniformly spaced reference
 * coordinates are generated and the rational approximant P(x)/Q(x) is
 * evaluated there.  The result is written as VTK_POLY_LINE cells so
 * Paraview renders a smooth curve.
 *
 * @param mesh          DG mesh
 * @param ref           Reference element
 * @param u             DG solution FieldVector (GLL nodal values)
 * @param solver        PadeLegendreSolver<Real,Index> constructed for the desired [L/M]
 * @param filename      Output VTK file path
 * @param fieldName     Scalar field name in the VTK output
 * @param t             Current simulation time
 * @param refineFactor  Points per original GLL node (≥ 1, default 4)
 */
void writePadeVTK(const Mesh<Real, Device, Index>&         mesh,
                  const ReferenceElement<Real, Index>&     ref,
                  const FieldVector<Real, Device, Index>&  u,
                  const PadeLegendreSolver<Real, Index>&   solver,
                  const std::string&                       filename,
                  const std::string&                       fieldName   = "u_pade",
                  Real                                     t           = Real(0),
                  Index                                    refineFactor = Index(4))
{
  if (refineFactor < 1)
    throw std::invalid_argument("writePadeVTK: refineFactor must be >= 1");

  const Index K        = mesh.numElements();
  const Index Np       = ref.numDOF();
  const Index nFine    = refineFactor * Np;   // fine points per element
  const Index totalPts = K * nFine;

  // Build equidistant reference coordinates once
  std::vector<Real> ref_pts(nFine);
  for (Index i = 0; i < nFine; ++i)
    ref_pts[i] = Real(-1) + Real(2) * i / Real(nFine - 1);

  std::ofstream f(filename);
  if (!f)
    throw std::runtime_error("writePadeVTK: cannot open '" + filename + "'");

  f << std::scientific << std::setprecision(10);

  // ── VTK header ────────────────────────────────────────────────────────────
  f << "# vtk DataFile Version 3.0\n";
  f << "DG1D Pade-Legendre reconstruction: " << fieldName << " t=" << t << "\n";
  f << "ASCII\n";
  f << "DATASET UNSTRUCTURED_GRID\n";

  // ── Points ────────────────────────────────────────────────────────────────
  f << "POINTS " << totalPts << " double\n";
  for (Index k = 0; k < K; ++k)
  {
    const Real xL = mesh.leftVertex(k);
    const Real h  = mesh.elementSize(k);
    for (Index i = 0; i < nFine; ++i)
    {
      const Real r = ref_pts[i];
      const Real x = xL + (r + Real(1)) * Real(0.5) * h;
      f << x << " 0.0 0.0\n";
    }
  }

  // ── Cells ─────────────────────────────────────────────────────────────────
  const Index cellListSize = K * (nFine + 1);
  f << "\nCELLS " << K << " " << cellListSize << "\n";
  for (Index k = 0; k < K; ++k)
  {
    f << nFine;
    for (Index i = 0; i < nFine; ++i)
      f << " " << k * nFine + i;
    f << "\n";
  }

  f << "\nCELL_TYPES " << K << "\n";
  for (Index k = 0; k < K; ++k)
    f << "4\n";

  // ── Point data ────────────────────────────────────────────────────────────
  f << "\nPOINT_DATA " << totalPts << "\n";
  f << "SCALARS " << fieldName << " double 1\n";
  f << "LOOKUP_TABLE default\n";
  for (Index k = 0; k < K; ++k)
  {
    const Real* uk  = u.elementPtr(k);
    const auto  vals = solver.reconstruct(uk, ref_pts);
    for (Index i = 0; i < nFine; ++i)
      f << vals[i] << "\n";
  }
}

// ─────────────────────────────────────────────────────────────────────────────
//  writePadeTimeSeriesVTK
//
//  Convenience wrapper: produces numbered VTK frames for the Padé reconstruction
//  on the finer equidistant grid, suitable for loading as a time series in
//  Paraview alongside the regular output frames.
//
//  Produces: basename_000000.vtk, basename_000001.vtk, ...
// ─────────────────────────────────────────────────────────────────────────────
template
<
  class Real,
  class Device,
  class Index
>
/**
 * @brief Convenience wrapper to write numbered Padé VTK frames for a time series.
 *
 * Produces files named basename_XXXXXX.vtk (6-digit zero-padded) that can be
 * loaded as a grouped time series in Paraview and compared directly with the
 * standard DG output series.
 */
void writePadeTimeSeriesVTK(const Mesh<Real, Device, Index>&        mesh,
                             const ReferenceElement<Real, Index>&    ref,
                             const FieldVector<Real, Device, Index>& u,
                             const PadeLegendreSolver<Real, Index>&  solver,
                             const std::string&                      basename,
                             int                                     frameIndex,
                             Real                                    t,
                             const std::string&                      fieldName    = "u_pade",
                             Index                                   refineFactor = Index(4))
{
  std::ostringstream oss;
  oss << basename << "_" << std::setfill('0') << std::setw(6) << frameIndex << ".vtk";
  writePadeVTK(mesh, ref, u, solver, oss.str(), fieldName, t, refineFactor);
}

} // namespace DG
