#pragma once

#include <print>

#include "FieldVector.hpp"
#include "Integrator.hpp"
#include "Operator.hpp"
#include "ReferenceElement.hpp"
#include "Traits.h"

// types defined in Traits
using RealType     = typename Traits::RealType;
using Device       = typename Traits::Device;
using IndexType    = typename Traits::IndexType;
using MeshType     = typename Traits::MeshType;
using MyMesh       = typename Traits::MyMesh;
using OperatorType = typename Traits::Operator;

namespace TNL::DGM {

/**
 * @class Solver
 * @brief Manages all the parts of the DGM solver
 * @template Model - physical model
 */
template < class Model >
class Solver
{
public:

   Solver() = default;

   static void init(int argc, char* argv[])
   {
      // 1. read parameter container
      // 2. set parameters
      // 3. read and load mesh
      // 4. initialize mesh data
      // 5. construct ReferenceElement
      // 6. impose initial conditions
   }

   /**
    * @brief Initialize the mesh
    */
   void initMesh() {}

   /**
    * @brief Helper to initialize the mesh
    */
   void initMeshData() {}

   /**
    * @brief Compute time step for ODE solver
    * TODO: switch to TNL::ODESolver and its time step computation
    */
   void compute_time_step()
   {
      // parallel reduction
      RealType r_min = 2.0;
      RealType max_speed = 0.0;
      for (int k = 0; k < Mesh.numElements(); k++)
      {
        for (int i = 1; i < RefElement.numDOF(); i++)
        {
          r_min = TNL::min(r_min, RefElement.nodes()[i] - RefElement.nodes()[i-1]);
          max_speed = TNL::argAbsMax(max_speed, advection_speed(Data.elementPtr(k)[i]));
        }
      }
      RealType x_min = r_min * Mesh.minJacobian();

      time_step_ = ODESolver.computeDt(x_min, max_speed, PolynomialOrder);
   }

   void setTime(RealType time)
   {
      time_ = time;
   }

   void updateTime(RealType time_increment)
   {
      time_ += time_increment;
   }

   /**
    * @brief Integrate the ODE system from t_0 to T
    */
   void integrate()
   {

   }

   /**
    * @brief Advance one step in time
    * TODO: Switch to TNL::ODESolver
    */
   void step()
   {
      ODESolver.step(Data, time_step_, time_);
      updateTime(time_step_);
   }

   /**
    * @brief Apply filter to the data
    * TODO: Implement
    */
   void applyFilter() {}

   /**
    * @brief Make a snapshot of the solution
    */
   void makeSnapshot()
   {
      DG::writeTimeSeriesVTK(Mesh, RefElement, Data, "output/output", frame++, time_);
   }

   /**
    * @brief Helper that reads the user config files
    */
   void readUserConfig() {}

   /**
    * @brief Logging function
    */
   void writeLog()
   {
      std::println("Starting simulation with: ");
      std::println("\tK = {}", Mesh.numElements());
      std::println("\tN = {}", RefElement.numDOF());
      std::println("\tdt = {}", time_step_);
   }

   /**
    * @brief Impose initial conditions
    * TODO: implement
    */
   template < class lambda_function >
   void imposeIC(lambda_function initial_conditions) {}

   FieldVector< RealType, Device, IndexType >& getData() { return Data; }

   MyMesh& getMesh() { return Mesh; }

   ReferenceElement< RealType, IndexType >& getReferenceElement() { return RefElement; }

private:

   // this should be modifiable from config file? or Traits?
   SSPRK< RealType, Device >                  ODESolver;
   // will become mesh data
   FieldVector< RealType, Device, IndexType > Data;
   // this should become MeshType
   MyMesh                                     Mesh;

   ReferenceElement< RealType, IndexType >    RefElement;
   OperatorType                               RHSOperator;

   // TODO
   DG::RoeFlux<RealType> numerical_flux( advection_speed, physical_flux );

   RealType                                   time_;
   RealType                                   time_step_;
   // used when taking a snapshot
   RealType                                   frame;
};

} // TNL::DGM
