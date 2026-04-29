#pragma once

#include <TNL/Timer.h>
#include <iostream>
#include "Traits.h"

#include <TNL/Containers/Expressions/ExpressionTemplates.h>

// ------------------------------------------------------------------------------------------------------------------ //

template< typename Traits >
class Simulation
{
   using GlobalIndexType = typename Traits::GlobalIndexType;
   using LocalIndexType = typename Traits::LocalIndexType;
   using IndexArrayType = typename Traits::IndexArrayType;
   using IntegerType = typename Traits::IntegerType;
   using VectorType = typename Traits::VectorType;
   using WriterType = typename Traits::WriterType;
   using MeshType = typename Traits::MeshType;
   using RealType = typename Traits::RealType;
   using Device = typename Traits::Device;

private:
   // boundary conditions
   RealType rho_0;
   RealType p_0;
   RealType alpha_0;
   RealType p_outlet;
   RealType u_0;
   RealType v_0;

   // initial conditions

public:
   // attributes
   const MeshData< Traits > meshData;

   // --------------------------------------------------------------------------------------------------------------- //

   // constructor
   Simulation( const std::string& fileName )
   : meshData( fileName )
   {
      printf( "> Simulation constructed.\n" );
   }

   // --------------------------------------------------------------------------------------------------------------- //

   // BoundaryConditions bcs;

   // InitialConditions ics;

   // methods
   void
   run( const RealType CFL, const RealType epsilon, const IntegerType repsMax )
   {
      TNL::Timer t_computeLocalTimeSteps;
      TNL::Timer t_updateBoundaries;
      TNL::Timer t_updatePoints;
      TNL::Timer t_reconstructValues;
      TNL::Timer t_computeFluxes;
      TNL::Timer t_computeLambdaStars;
      TNL::Timer t_computeDiagonalCoeffs;
      TNL::Timer t_forwardSweep;
      TNL::Timer t_backwardSweep;
      TNL::Timer t_computeRezi;
      RealType t_totalTime = 0;


      RealType rezi = 1;
      IntegerType reps = 0;
      while( rezi > epsilon && reps < repsMax ) {
         t_computeLocalTimeSteps.start();
         computeLocalTimeSteps( CFL );
         setGlobalTimeStep();
         t_computeLocalTimeSteps.stop();

         t_updateBoundaries.start();
         updateBoundaries();
         t_updateBoundaries.stop();

         t_updatePoints.start();
         updatePoints();
         t_updatePoints.stop();

         t_reconstructValues.start();
         reconstructValues();
         t_reconstructValues.stop();

         t_computeFluxes.start();
         computeFluxes();
         t_computeFluxes.stop();

         // LU-SGS
         t_computeLambdaStars.start();
         computeLambdaStars();
         t_computeLambdaStars.stop();

         t_computeDiagonalCoeffs.start();
         computeDiagonalCoeffs();
         t_computeDiagonalCoeffs.stop();

         t_forwardSweep.start();
         forwardSweep();
         t_forwardSweep.stop();

         t_backwardSweep.start();
         backwardSweep();
         t_backwardSweep.stop();
         // \LU-SGS

         t_computeRezi.start();
         rezi = computeRezi();
         t_computeRezi.stop();

         fluidData.w_old = fluidData.w_new;
         reps++;
      }

      t_totalTime = t_computeLocalTimeSteps.getRealTime() + t_updateBoundaries.getRealTime() + t_reconstructValues.getRealTime()
                  + t_updatePoints.getRealTime() + t_computeFluxes.getRealTime() + t_computeLambdaStars.getRealTime()
                  + t_computeDiagonalCoeffs.getRealTime() + t_forwardSweep.getRealTime() + t_backwardSweep.getRealTime()
                  + t_computeRezi.getRealTime();

      printf( " ** rezi = %f ** \n", rezi );

      std::cout << std::fixed << std::setprecision( 3 );
      std::cout << "totalTime......................: (" << std::setw( 7 ) << t_totalTime / t_totalTime * 100 << " %) "
                << t_totalTime << " sec." << std::endl;

      std::cout << "-------------------------------------------------------" << std::endl;

      std::cout << "t_computeLocalTimeSteps........: (" << std::setw( 7 )
                << t_computeLocalTimeSteps.getRealTime() / t_totalTime * 100 << " %) " << t_computeLocalTimeSteps.getRealTime()
                << " sec." << std::endl;

      std::cout << "t_updateBoundaries.............: (" << std::setw( 7 )
                << t_updateBoundaries.getRealTime() / t_totalTime * 100 << " %) " << t_updateBoundaries.getRealTime() << " sec."
                << std::endl;

      std::cout << "t_updatePoints.................: (" << std::setw( 7 ) << t_updatePoints.getRealTime() / t_totalTime * 100
                << " %) " << t_updatePoints.getRealTime() << " sec." << std::endl;

      std::cout << "t_reconstructValues............: (" << std::setw( 7 )
                << t_reconstructValues.getRealTime() / t_totalTime * 100 << " %) " << t_reconstructValues.getRealTime()
                << " sec." << std::endl;

      std::cout << "t_computeFluxes................: (" << std::setw( 7 ) << t_computeFluxes.getRealTime() / t_totalTime * 100
                << " %) " << t_computeFluxes.getRealTime() << " sec." << std::endl;

      std::cout << "t_computeLambdaStars...........: (" << std::setw( 7 )
                << t_computeLambdaStars.getRealTime() / t_totalTime * 100 << " %) " << t_computeLambdaStars.getRealTime()
                << " sec." << std::endl;

      std::cout << "t_computeDiagonalCoeffs........: (" << std::setw( 7 )
                << t_computeDiagonalCoeffs.getRealTime() / t_totalTime * 100 << " %) " << t_computeDiagonalCoeffs.getRealTime()
                << " sec." << std::endl;

      std::cout << "t_forwardSweep.................: (" << std::setw( 7 ) << t_forwardSweep.getRealTime() / t_totalTime * 100
                << " %) " << t_forwardSweep.getRealTime() << " sec." << std::endl;

      std::cout << "t_backwardSweep................: (" << std::setw( 7 ) << t_backwardSweep.getRealTime() / t_totalTime * 100
                << " %) " << t_backwardSweep.getRealTime() << " sec." << std::endl;

      std::cout << "t_computeRezi..................: (" << std::setw( 7 ) << t_computeRezi.getRealTime() / t_totalTime * 100
                << " %) " << t_computeRezi.getRealTime() << " sec." << std::endl;
   }

   // --------------------------------------------------------------------------------------------------------------- //

   void
   setTestInitialCondition()
   {
      auto w_old_view = fluidData.w_old.getView();

      auto setInitialCondition = [ = ] __cuda_callable__( const GlobalIndexType iterIndex ) mutable
      {
         const auto DIMENSION = MeshType::getMeshDimension();
         const auto cell = meshData.meshPointer->template getEntity< DIMENSION >( iterIndex );
         const auto cellCenter = getEntityCenter( *meshData.meshPointer, cell );
         ConservativeType wInitial{};
         wInitial[ 0 ] = 1;
         wInitial[ 1 ] = cellCenter.x();
         wInitial[ 2 ] = cellCenter.y();
         wInitial[ 3 ] = 2;
         w_old_view[ iterIndex ] = wInitial;
         printf( "> processed cell %d at [%f, %f]\n", iterIndex, cellCenter.x(), cellCenter.y() );
      };
      fluidData.w_new = fluidData.w_old;
      fluidData.w_one = fluidData.w_old;
      parallelFor< Device >( 0, w_old_view.getSize(), setInitialCondition );
   }

   // --------------------------------------------------------------------------------------------------------------- //

   void
   setInitialCondition( const ConservativeType& wInitial )
   {
      auto w_old_view = fluidData.w_old.getView();
      auto w_one_view = fluidData.w_one.getView();
      auto w_new_view = fluidData.w_new.getView();

      auto setInitialCondition = [ = ] __cuda_callable__( const GlobalIndexType iterIndex ) mutable
      {
         w_old_view[ iterIndex ] = wInitial;
         w_one_view[ iterIndex ] = wInitial;
         w_new_view[ iterIndex ] = wInitial;
      };
      parallelFor< Device >( 0, w_old_view.getSize(), setInitialCondition );
   }

   // --------------------------------------------------------------------------------------------------------------- //

   void
   computeLocalTimeSteps( const RealType CFL )
   {
      const auto innerCellDTVectors_const_view = meshData.primary.innerCellDTVectors.getConstView();
      const auto innerCellIDs_const_view = meshData.IDs.innerCellIDs.getConstView();
      const auto w_old_const_view = fluidData.w_old.getConstView();
      auto innerCellDTs_view = fluidData.innerCellDTs.getView();

      auto computeSingleTimeStep = [ = ] __cuda_callable__( const GlobalIndexType iterIndex ) mutable
      {
         const GlobalIndexType cellIndex = innerCellIDs_const_view[ iterIndex ];

         const Primitive< Traits > currPV( w_old_const_view[ cellIndex ] );
         const RealType u = currPV.u;
         const RealType v = currPV.v;
         const RealType c = currPV.c;

         const VectorType xi = innerCellDTVectors_const_view( iterIndex, 0 );
         const VectorType eta = innerCellDTVectors_const_view( iterIndex, 1 );

         const RealType xi_norm = TNL::l2Norm( xi );
         const RealType eta_norm = TNL::l2Norm( eta );

         const RealType u_xi = fabs( ( u * xi.x() + v * xi.y() ) / xi_norm );
         const RealType u_eta = fabs( ( u * eta.x() + v * eta.y() ) / eta_norm );

         const RealType d_xi = ( u_xi + c ) / xi_norm;
         const RealType d_eta = ( u_eta + c ) / eta_norm;

         const RealType res = CFL / ( d_xi + d_eta );
         innerCellDTs_view[ iterIndex ] = res;
      };
      parallelFor< Device >( 0, meshData.IDs.innerCellIDs.getSize(), computeSingleTimeStep );
   }

   // --------------------------------------------------------------------------------------------------------------- //

   void
   setGlobalTimeStep()
   {
      auto innerCellDTs_view = fluidData.innerCellDTs.getView();

      const RealType maxTimeStep = TNL::max( fluidData.innerCellDTs );
      auto assignTimeStep = [ = ] __cuda_callable__( const GlobalIndexType iterIndex ) mutable
      {
         innerCellDTs_view[ iterIndex ] = maxTimeStep;
      };
      parallelFor< Device >( 0, meshData.IDs.innerCellIDs.getSize(), assignTimeStep );
   }

   // --------------------------------------------------------------------------------------------------------------- //

   void
   updatePoints()
   {
      const auto pointNeighborOffsets_const_view = meshData.connectivity.pointNeighborOffsets.getConstView();
      const auto pointNeighborCounts_const_view = meshData.connectivity.pointNeighborCounts.getConstView();
      const auto pointNeighborIDs_const_view = meshData.connectivity.pointNeighborIDs.getConstView();
      const auto w_old_const_view = fluidData.w_old.getConstView();
      auto w_points_view = fluidData.w_points.getView();

      auto processPoint = [ = ] __cuda_callable__( const GlobalIndexType iterIndex ) mutable
      {
         ConservativeType w_currPoint{};
         RealType weightSum = 0;
         const IntegerType offset = pointNeighborOffsets_const_view[ iterIndex ];
         const IntegerType neighborCount = pointNeighborCounts_const_view[ iterIndex ];
         for( IntegerType neighborIndex = 0; neighborIndex < neighborCount; ++neighborIndex ) {
            const auto cellNeighborIndex = pointNeighborIDs_const_view[ offset + neighborIndex ];
            const auto wCell = w_old_const_view[ cellNeighborIndex ];
            w_currPoint += wCell;
         }
         w_currPoint /= neighborCount;
         w_points_view[ iterIndex ] = w_currPoint;
      };
      parallelFor< Device >( 0, w_points_view.getSize(), processPoint );
   }

   // --------------------------------------------------------------------------------------------------------------- //

   void
   reconstructValues()
   {
      // insert modular function
      // first order - trivial reconstruction
      reconstruction.reconstructValues( fluidData.w_old, fluidData.wi_reconstructed, fluidData.wj_reconstructed );
   }
   // --------------------------------------------------------------------------------------------------------------- //

   void
   computeFluxes()
   {
      const auto innerCellInterfacePoints_const_view = meshData.connectivity.innerCellInterfacePoints.getConstView();
      const auto innerCellNeighborOffsets_const_view = meshData.connectivity.innerCellNeighborOffsets.getConstView();
      const auto innerCellNeighborCounts_const_view = meshData.connectivity.innerCellNeighborCounts.getConstView();
      const auto innerCellNeighborIDs_const_view = meshData.connectivity.innerCellNeighborIDs.getConstView();
      const auto cellNeighborOffsets_const_view = meshData.connectivity.cellNeighborOffsets.getConstView();
      const auto cellNeighborCounts_const_view = meshData.connectivity.cellNeighborCounts.getConstView();
      const auto cellNeighborIDs_const_view = meshData.connectivity.cellNeighborIDs.getConstView();

      const auto innerCellScaledNormals_const_view = meshData.primary.innerCellScaledNormals.getConstView();
      const auto dualMeshScaledNormals_const_view = meshData.dual.dualMeshScaledNormals.getConstView();
      const auto dualMeshVolumes_const_view = meshData.dual.dualMeshVolumes.getConstView();
      const auto innerCellIDs_const_view = meshData.IDs.innerCellIDs.getConstView();

      const auto wi_reconstructed_const_view = fluidData.wi_reconstructed.getConstView();
      const auto wj_reconstructed_const_view = fluidData.wj_reconstructed.getConstView();
      const auto w_points_const_view = fluidData.w_points.getConstView();
      const auto w_old_const_view = fluidData.w_old.getConstView();
      auto innerCellR_view = fluidData.innerCellR.getView();

      auto computeFlux = [ = ] __cuda_callable__( const GlobalIndexType iterIndex ) mutable
      {
         const auto cellIndex = innerCellIDs_const_view[ iterIndex ];

         ConservativeType convectiveTermsSum{};
         ConservativeType diffusiveTermsSum{};
         const auto innerOffset = innerCellNeighborOffsets_const_view[ iterIndex ];
         const auto globalOffset = cellNeighborOffsets_const_view[ cellIndex ];
         const auto neighborCount = innerCellNeighborCounts_const_view[ iterIndex ];
         const auto wi = w_old_const_view[ cellIndex ];
         for( LocalIndexType neighborIndex = 0; neighborIndex < neighborCount; ++neighborIndex ) {
            const auto scaledNormal = innerCellScaledNormals_const_view[ innerOffset + neighborIndex ];
            const auto wi_reconstructed = wi_reconstructed_const_view[ globalOffset + neighborIndex ];
            const auto wj_reconstructed = wj_reconstructed_const_view[ globalOffset + neighborIndex ];
            const auto convectiveIncrement = HLLC< Traits >( wi_reconstructed, wj_reconstructed, scaledNormal );
            convectiveTermsSum += convectiveIncrement;

            const auto cellNeighborIndex = innerCellNeighborIDs_const_view[ innerOffset + neighborIndex ];
            const auto vertexIndexR = innerCellInterfacePoints_const_view( innerOffset + neighborIndex, 0 );
            const auto vertexIndexL = innerCellInterfacePoints_const_view( innerOffset + neighborIndex, 1 );
            const auto wj = w_old_const_view[ cellNeighborIndex ];
            const auto wr = w_points_const_view[ vertexIndexR ];
            const auto wl = w_points_const_view[ vertexIndexL ];
            const auto sn1 = dualMeshScaledNormals_const_view( innerOffset + neighborIndex, 0 );
            const auto sn2 = dualMeshScaledNormals_const_view( innerOffset + neighborIndex, 1 );
            const auto sn3 = dualMeshScaledNormals_const_view( innerOffset + neighborIndex, 2 );
            const auto sn4 = dualMeshScaledNormals_const_view( innerOffset + neighborIndex, 3 );
            const auto dualMeshVolume = dualMeshVolumes_const_view[ innerOffset + neighborIndex ];
            const auto diffusiveIncrement =
               computeDiffusiveTerms( wi, wj, wr, wl, scaledNormal, sn1, sn2, sn3, sn4, dualMeshVolume );
            diffusiveTermsSum += diffusiveIncrement;
         }
         innerCellR_view[ iterIndex ] = convectiveTermsSum - diffusiveTermsSum;
      };
      parallelFor< Device >( 0, meshData.IDs.innerCellIDs.getSize(), computeFlux );
   }

   // --------------------------------------------------------------------------------------------------------------- //

   void
   updateDomainExplicit()
   {
      const auto innerCellAreas_const_view = meshData.primary.innerCellAreas.getConstView();
      const auto innerCellDTs_const_view = fluidData.innerCellDTs.getConstView();
      const auto innerCellIDs_const_view = meshData.IDs.innerCellIDs.getConstView();
      const auto w_old_const_view = fluidData.w_old.getConstView();
      auto innerCellR_view = fluidData.innerCellR.getView();
      auto w_new_view = fluidData.w_new.getView();

      auto updateCell = [ = ] __cuda_callable__( const GlobalIndexType iterIndex ) mutable
      {
         const RealType cellArea = innerCellAreas_const_view[ iterIndex ];
         const RealType cellDT = innerCellDTs_const_view[ iterIndex ];
         const GlobalIndexType cellIndex = innerCellIDs_const_view[ iterIndex ];
         w_new_view[ cellIndex ] = w_old_const_view[ cellIndex ] - innerCellR_view[ iterIndex ] * cellDT / cellArea;
      };
      parallelFor< Device >( 0, meshData.IDs.innerCellIDs.getSize(), updateCell );
   }

   // --------------------------------------------------------------------------------------------------------------- //

   RealType
   computeRezi()
   {
      const auto innerCellIDs_const_view = meshData.IDs.innerCellIDs.getConstView();
      const auto innerCellAreas_const_view = meshData.primary.innerCellAreas.getConstView();
      const auto innerCellDTs_const_view = fluidData.innerCellDTs.getConstView();
      const auto w_old_const_view = fluidData.w_old.getConstView();
      const auto w_new_const_view = fluidData.w_new.getConstView();

      // get value from one cell
      auto fetch = [ = ] __cuda_callable__( const GlobalIndexType iterIndex )
      {
         const GlobalIndexType cellIndex = innerCellIDs_const_view[ iterIndex ];
         const RealType cellArea = innerCellAreas_const_view[ iterIndex ];
         const RealType cellDT = innerCellDTs_const_view[ iterIndex ];
         const RealType rho_new = w_new_const_view[ cellIndex ].getElement( 0 );
         const RealType rho_old = w_old_const_view[ cellIndex ].getElement( 0 );

         return cellArea * pow( ( rho_new - rho_old ) / cellDT, 2 );
      };

      auto reduction = [ = ] __cuda_callable__( const RealType& sum, const RealType& increment )
      {
         return sum + increment;
      };

      // parallel reduction
      return log( sqrt( reduce< Device >( 0, meshData.IDs.innerCellIDs.getSize(), fetch, reduction, 0.0 ) ) );
   }

   // --------------------------------------------------------------------------------------------------------------- //

   void
   exportResults( const std::string& outputFileName )
   {
      Vector< RealType, Device > res_mach( fluidData.w_points.getSize() );
      Vector< RealType, Device > res_cp( fluidData.w_points.getSize() );
      Vector< RealType, Device > res_u( fluidData.w_points.getSize() );
      Vector< RealType, Device > res_v( fluidData.w_points.getSize() );
      Vector< RealType, Device > res_p( fluidData.w_points.getSize() );

      const auto w_points_const_view = fluidData.w_points.getConstView();
      auto res_mach_view = res_mach.getView();
      auto res_cp_view = res_cp.getView();
      auto res_u_view = res_u.getView();
      auto res_v_view = res_v.getView();
      auto res_p_view = res_p.getView();

      const RealType rho_0_device = rho_0;
      const RealType p_0_device = p_0;

      auto writeData = [ = ] __cuda_callable__( const GlobalIndexType iterIndex ) mutable
      {
         const Primitive< Traits > pv_ref( w_points_const_view[ 0 ] );
         const Primitive< Traits > pv( w_points_const_view[ iterIndex ] );
         res_mach_view[ iterIndex ] = getMachNumber( pv.U, pv.c );
         res_cp_view[ iterIndex ] = getCP( pv_ref.p, pv.p, rho_0_device, p_0_device );
         res_u_view[ iterIndex ] = pv.u;
         res_v_view[ iterIndex ] = pv.v;
         res_p_view[ iterIndex ] = pv.p;
      };
      parallelFor< Device >( 0, meshData.connectivity.pointNeighborCounts.getSize(), writeData );

      IndexArrayType iterIDs( meshData.cellEntityIDs.getSize() );
      auto iterIDs_view = iterIDs.getView();
      auto writeID = [ = ] __cuda_callable__( const GlobalIndexType iterIndex ) mutable
      {
         iterIDs_view[ iterIndex ] = iterIndex;
      };
      parallelFor< Device >( 0, iterIDs.getSize(), writeID );

      // bring data over to host
      Vector< IntegerType, TNL::Devices::Host > cellEntityIDs_host( meshData.cellEntityIDs );
      Vector< RealType, TNL::Devices::Host > res_mach_host( res_mach );
      Vector< RealType, TNL::Devices::Host > res_cp_host( res_cp );
      Vector< RealType, TNL::Devices::Host > res_u_host( res_u );
      Vector< RealType, TNL::Devices::Host > res_v_host( res_v );
      Vector< RealType, TNL::Devices::Host > res_p_host( res_p );

      Vector< GlobalIndexType, TNL::Devices::Host > iterIDs_host( iterIDs );

      std::ofstream fileToSaveMesh( outputFileName );
      WriterType writer( fileToSaveMesh );
      writer.writeEntities( meshData.hostMesh );
      writer.template writeCellData< Vector< IntegerType, TNL::Devices::Host > >( iterIDs_host, "iterIDs", 1 );
      writer.template writeCellData< Vector< IntegerType, TNL::Devices::Host > >( cellEntityIDs_host, "cell-entity-IDs", 1 );
      writer.template writePointData< Vector< RealType, TNL::Devices::Host > >( res_cp_host, "pressure-coefficient", 1 );
      writer.template writePointData< Vector< RealType, TNL::Devices::Host > >( res_mach_host, "mach-number", 1 );
      writer.template writePointData< Vector< RealType, TNL::Devices::Host > >( res_u_host, "u", 1 );
      writer.template writePointData< Vector< RealType, TNL::Devices::Host > >( res_v_host, "v", 1 );
      writer.template writePointData< Vector< RealType, TNL::Devices::Host > >( res_p_host, "pressure", 1 );
      printf( "Data written successfully.\n" );
   }

   // --------------------------------------------------------------------------------------------------------------- //
};

// ------------------------------------------------------------------------------------------------------------------ //
