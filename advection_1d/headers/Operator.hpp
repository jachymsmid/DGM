#pragma once

#include "Mesh.hpp"
#include "ReferenceElement.hpp"
#include "FieldVector.hpp"
#include "NumericalFlux.hpp"

namespace DG {

template
<
  class Real = double,
  class Device = TNL::Devices::Host,
  class Index = int
>
class Operator
{
public:
  using MeshType = Mesh<Real, Device, Index>;
  using Ref  = ReferenceElement<Real, Index>;
  using Field = FieldVector<Real, Device, Index>;

  // Constructor
  Operator(
      const MeshType& mesh,
      const Ref& ref,
      const NumericalFlux<Real>& flux,
      std::function<Real(Real)> physFlux)
      : mesh_(mesh), ref_(ref), flux_(flux), physFlux_(physFlux)
  {}

  // return the rhs operatro as a std::function
  std::function<void(const Field&, Field&, const Real& time)> rhsFunction() const
  {
    return [this](const Field& u, Field& rhs, const Real& time) { computeRHS(u, rhs, time); };
  }

  // compute rhs = du/dt into 'res' given current state 'u' and 'time'
  void computeRHS(const Field& u, Field& res, const Real& time) const
  {
    const Index K  = mesh_.numElements();
    const Index Np = ref_.numDOF();
    const auto& Dr   = ref_.Dr();
    const auto& LIFT = ref_.LIFT();


    std::vector<Real> vol(Np), flocal(Np), fluxJump(2);

    // for every mesh element
    for (Index k = 0; k < K; ++k)
    {
      // create cell arrays
      const Real* uk = u.elementPtr(k);
      Real* rk = res.elementPtr(k);

      // physical flux at all nodes of the element
      for (Index i = 0; i < Np; ++i)
      {
        flocal[i] = physFlux_(uk[i]);
      }

      // --- volume term ---
      // vol = Dr * flocal
      for (Index i = 0; i < Np; ++i)
      {
        Real s = 0;
        for (Index j = 0; j < Np; ++j)
        {
          s += Dr(i,j) * flocal[j];
        }
        vol[i] = s;
      }

      // --- surface term ---
      // left face
      {
        // current element's value at left face
        Real u_int = physFlux_(uk[0]);
        // neighbour's face value
        Real u_ext;
        // BC
        if (mesh_.isBoundaryFace(k))
        {
          // periodic
          u_ext = physFlux_(u.elementPtr(K-1)[Np-1]);
          // ghost cell
          // u_ext = u_int;
        }
        else
        {
          u_ext = physFlux_(u.elementPtr(mesh_.leftCellOfFace(k))[Np-1]);
        }
        // fluxJump[0]  = fStar
        fluxJump[0]  = flux_.compute(u_ext, u_int, mesh_.leftNormal());
      }

      // right face
      {
        // current element's value at right face
        Real u_int = physFlux_(uk[Np-1]);
        // neighbour's face value
        Real u_ext;

        // BC
        if (mesh_.isBoundaryFace(k+1))
        {
          // periodic: right neighbour is first element's left end
          u_ext = physFlux_(u.elementPtr(0)[0]);
          // ghost cell: right neighbour takes its own value
          // u_ext = u_int;
        }
        else
        {
          u_ext = physFlux_(u.elementPtr(mesh_.rightCellOfFace(k+1))[0]);
        }
        // numerical flux f*(u-, u+)
        fluxJump[1]  = flux_.compute(u_int, u_ext, mesh_.rightNormal());
      }

      // --- assemble rhs ---
      // rhs = J^{-1} * (-vol + LIFT * fluxJump)
      Real Jinv = Real(1) / mesh_.jacobian(k);
      for (Index i = 0; i < Np; ++i)
      {
        Real lift = LIFT(i,0)*fluxJump[0] + LIFT(i,1)*fluxJump[1];
        rk[i] = Jinv * ( vol[i] + lift);
      }
    }
  }

private:
    const MeshType& mesh_;
    const Ref& ref_;
    const NumericalFlux<Real>& flux_;
    std::function<Real(Real)> physFlux_;
};

} // namespace DG
