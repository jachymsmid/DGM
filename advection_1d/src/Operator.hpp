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

  // compute rhs = du/dt into 'res' given current state 'u'.
  void computeRHS(const Field& u, Field& res) const
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
      for (Index i = 0; i < Np; ++i) flocal[i] = physFlux_(uk[i]);

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
        Real u_int = uk[0]; 
        // neighbour's face value
        Real u_ext; 
        // periodic BC
        // TODO: implement different BC
        if (mesh_.isBoundaryFace(k))
        {
          u_ext = u.elementPtr(K-1)[Np-1];
        }
        else
        {
          u_ext = u.elementPtr(mesh_.leftCellOfFace(k))[Np-1];
        }
        // numerical flux f*
        Real fStar = flux_.compute(u_ext, u_int, mesh_.leftNormal());
        // f(u_int) = f(u+)
        Real fLocal = physFlux_(u_int);
        fluxJump[0]  = mesh_.leftNormal() * (fStar - fLocal);
        // fluxJump[0]  = (fStar - fLocal);
      }

      // right face
      {
        // current element's value at right face
        Real u_int  = uk[Np-1];
        // neighbour's face value
        Real u_ext;
        // periodic BC
        if (mesh_.isBoundaryFace(k+1))
        {
          // periodic: right neighbour is first element's left end
          u_ext = u.elementPtr(0)[0];
        }
        else
        {
          u_ext = u.elementPtr(mesh_.rightCellOfFace(k+1))[0];
        }
        // numerical flux f*(u-, u+)
        Real fStar   = flux_.compute(u_int, u_ext, mesh_.rightNormal());
        // f(u_int) = f(u-)
        Real fLocal  = physFlux_(u_int);
        // contribution to rhs: n * (f* - f_local), n=+1 at right face
        fluxJump[1]  = mesh_.rightNormal() * (fStar - fLocal);
        // fluxJump[1]  = (fStar - fLocal);
      }

      // --- assemble rhs ---
      // rhs = J^{-1} * (-vol + LIFT * fluxJump)
      Real Jinv = Real(1) / mesh_.jacobian(k);
      for (Index i = 0; i < Np; ++i)
      {
        Real lift = LIFT(i,0)*fluxJump[0] + LIFT(i,1)*fluxJump[1];
        rk[i] = Jinv * (-vol[i] + lift);
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
