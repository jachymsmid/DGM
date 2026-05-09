#pragma once

#include "FieldVector.hpp"
#include "Integrator.hpp"
#include "ReferenceElement.hpp"
#include "Traits.h"

// types defined in Traits
using RealType  = typename Traits::RealType;
using Device    = typename Traits::Device;
using IndexType = typename Traits::IndexType;
using MeshType  = typename Traits::MeshType;

namespace TNL::DGM {

/**
 * @class Solver
 * @brief Manages all the parts of the DGM solver
 */
template < class Model >
class Solver
{
public:

   Solver() = default;

   void init(int argc, char* argv[])
   {
      // 1. read parameter container
      // 2. set parameters
      // 3. read and load mesh
      // 4. initialize mesh data
      // 5. impose initial conditions

   }

   void initMesh() {}

   void initMeshData() {}

   void compute_time_step()
   {
      // make a function for this, parallel reduction
      RealType r_min = 2.0;
      RealType max_speed = 0.0;
      for (int k = 0; k < mesh.numElements(); k++)
      {
        for (int i = 1; i < RefElement.numDOF(); i++)
        {
          r_min = TNL::min(r_min, RefElement.nodes()[i] - RefElement.nodes()[i-1]);
          max_speed = TNL::argAbsMax(max_speed, advection_speed(Data.elementPtr(k)[i]));
        }
      }
      RealType x_min = r_min * mesh.minJacobian();

      delta_t = ODESolver.computeDt(x_min, max_speed, PolynomialOrder);
   }

   void integrate()
   {

   }

   void makeSnapshot() {}

   void readUserConfig() {}

   void writeLog() {}

private:

   // this should be modifiable from config file? or Traits?
   SSPRK< RealType, Device >                  &ODESolver;
   // will become mesh data
   FieldVector< RealType, Device, IndexType > &Data;
   MeshType                                   &Mesh;
   RealType                                   delta_t;
   ReferenceElement<>                         &RefElement;
   IndexType                                  PolynomialOrder;
};

} // TNL::DGM

