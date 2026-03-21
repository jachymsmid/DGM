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
  using Mesh = Mesh<Real, Device, Index>;
  using Ref  = ReferenceElement<Real, Index>;
  using Field = FieldVector<Real, Device, Index>;

  // Constructor
  Operator(const Mesh& mesh, const Ref& ref,
             const NumericalFlux<Real>& flux,
             std::function<Real(Real)> physFlux)
      : mesh_(mesh), ref_(ref), flux_(flux), physFlux_(physFlux)
  {}

  // Compute rhs = du/dt into 'res' given current state 'u'.
  // Both must be allocated with (K, Np).
  void computeRHS(const Field& u, Field& res) const
  {
    const Index K  = mesh_.numElements();
    const Index Np = ref_.numDOF();
    const auto& Dr   = ref_.Dr();
    const auto& LIFT = ref_.LIFT();

    // --- 1. Volume term: vol_k = Dr * f(u_k)  ---
    std::vector<Real> vol(Np), flocal(Np), fluxJump(2);

    for (Index k = 0; k < K; ++k)
    {
      const Real* uk = u.elementPtr(k);
      Real* rk = res.elementPtr(k);

      // f at local DOFs
      for (Index i = 0; i < Np; ++i) flocal[i] = physFlux_(uk[i]);

      // vol = Dr * flocal
      for (Index i = 0; i < Np; ++i)
      {
        Real s = 0;
        for (Index j = 0; j < Np; ++j) s += Dr(i,j) * flocal[j];
        vol[i] = s;
      }

      // --- 2. Surface term: numerical flux at left (face k) and right (face k+1) ---
      // Left face (face index k, outward normal = -1)
     {
        Real u_int  = uk[0];          // current element's value at left face
        Real u_ext;                   // neighbour's value
        if (mesh_.isBoundaryFace(k)) {
            // periodic: left neighbour is last element's right end
            u_ext = u.elementPtr(K-1)[Np-1];
        } else {
            u_ext = u.elementPtr(mesh_.leftCellOfFace(k))[Np-1];
        }
        // upwind: a>0 so information comes from the left, u_minus=u_ext, u_plus=u_int
        Real fStar   = flux_.compute(u_ext, u_int, Real(-1));
        Real fLocal  = physFlux_(u_int);
        // contribution to rhs: n * (f* - f_local), n=-1 at left face
        fluxJump[0]  = Real(-1) * (fStar - fLocal);
      }

      // Right face (face index k+1, outward normal = +1)
      {
        Real u_int  = uk[Np-1];       // current element's value at right face
        Real u_ext;
        if (mesh_.isBoundaryFace(k+1)) {
            // periodic: right neighbour is first element's left end
            u_ext = u.elementPtr(0)[0];
        } else {
            u_ext = u.elementPtr(mesh_.rightCellOfFace(k+1))[0];
        }
        // upwind: u_minus=u_int (left), u_plus=u_ext (right)
        Real fStar   = flux_.compute(u_int, u_ext, Real(+1));
        Real fLocal  = physFlux_(u_int);
        // contribution to rhs: n * (f* - f_local), n=+1 at right face
        fluxJump[1]  = Real(+1) * (fStar - fLocal);
      }
      // --- 3. Assemble: rhs = J^{-1} * (-vol + LIFT * fluxJump) ---
      Real Jinv = Real(1) / mesh_.jacobian(k);
      for (Index i = 0; i < Np; ++i)
      {
        Real lift = LIFT(i,0)*fluxJump[0] + LIFT(i,1)*fluxJump[1];
        rk[i] = Jinv * (-vol[i] + lift);
      }
    }
  }

private:
    const Mesh& mesh_;
    const Ref& ref_;
    const NumericalFlux<Real>& flux_;
    std::function<Real(Real)> physFlux_;
};

} // namespace DG
