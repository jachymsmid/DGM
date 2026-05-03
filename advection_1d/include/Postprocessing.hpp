#pragma once

#include "FieldVector.hpp"
#include "ReferenceElement.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace TNL::DGM {

template< class Real = double, class Index = int >
struct PadeApproximant
{
  std::vector< Real > p_coeffs;
  std::vector< Real > q_coeffs;
  std::vector< Real > modal_coeffs;
  bool valid{ true };
  Index L{ 0 };
  Index M{ 0 };
};

template< class Real = double, class Index = int >
class PadeLegendreSolver
{
public:
  using ReferenceElementType = ReferenceElement< Real, Index >;
  using RefMatrix = typename ReferenceElementType::Matrix;
  using Field = FieldVector< Real, TNL::Devices::Host, Index >;
  using EigenMat = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >;
  using EigenVec = Eigen::Matrix< double, Eigen::Dynamic, 1 >;

  PadeLegendreSolver( const ReferenceElementType& ref,
                      Index L,
                      Index M,
                      Real fallback_tol = Real( 100 ) * std::numeric_limits< Real >::epsilon() )
  : ref_( ref ),
    L_( L ),
    M_( M ),
    N_( ref.order() ),
    Np_( ref.numDOF() ),
    fallback_tol_( fallback_tol ),
    vinv_( ref.Vinv() )
  {
    if( L_ < 0 || M_ < 0 )
      throw std::invalid_argument( "PadeLegendreSolver: L and M must be >= 0" );
    if( L_ + M_ > N_ )
      throw std::invalid_argument( "PadeLegendreSolver: L + M must be <= N" );
  }

  std::vector< Real > modalCoeffs( const Real* nodal_vals ) const
  {
    std::vector< Real > c( Np_, Real( 0 ) );
    for( Index n = 0; n < Np_; ++n ) {
      Real s = Real( 0 );
      for( Index i = 0; i < Np_; ++i )
        s += vinv_.getElement( n, i ) * nodal_vals[ i ];
      c[ n ] = s * std::sqrt( Real( 2 * n + 1 ) / Real( 2 ) );
    }
    return c;
  }

  PadeApproximant< Real, Index > buildApproximant( const std::vector< Real >& c ) const
  {
    if( static_cast< Index >( c.size() ) != Np_ )
      throw std::invalid_argument( "PadeLegendreSolver::buildApproximant: invalid modal size." );

    PadeApproximant< Real, Index > approx;
    approx.L = L_;
    approx.M = M_;
    approx.modal_coeffs = c;
    approx.p_coeffs.assign( L_ + 1, Real( 0 ) );
    approx.q_coeffs.assign( M_ + 1, Real( 0 ) );
    approx.q_coeffs[ 0 ] = Real( 1 );

    const EigenMat B = buildBMatrix_( c );

    if( M_ == 0 ) {
      for( Index k = 0; k <= L_; ++k )
        approx.p_coeffs[ k ] = static_cast< Real >( ( ( 2.0 * k + 1.0 ) / 2.0 ) * B( k, 0 ) );
      return approx;
    }

    const Index nRows = N_ - L_;
    EigenMat A( nRows, M_ );
    EigenVec b( nRows );

    for( Index k = L_ + 1; k <= N_; ++k ) {
      const Index row = k - L_ - 1;
      b( row ) = B( k, 0 );
      for( Index m = 1; m <= M_; ++m )
        A( row, m - 1 ) = B( k, m );
    }

    const EigenVec q_vec = A.colPivHouseholderQr().solve( -b );
    const double residual = ( A * q_vec + b ).norm();
    const double scale = std::max( 1.0, b.norm() );
    if( residual > 1e-10 * scale )
      approx.valid = false;

    for( Index m = 1; m <= M_; ++m )
      approx.q_coeffs[ m ] = static_cast< Real >( q_vec( m - 1 ) );

    for( Index k = 0; k <= L_; ++k ) {
      double s = 0.0;
      for( Index m = 0; m <= M_; ++m )
        s += static_cast< double >( approx.q_coeffs[ m ] ) * B( k, m );
      approx.p_coeffs[ k ] = static_cast< Real >( ( ( 2.0 * k + 1.0 ) / 2.0 ) * s );
    }

    return approx;
  }

  Real evaluate( const PadeApproximant< Real, Index >& approx, Real x ) const
  {
    if( ! approx.valid )
      return polynomialFromModal_( approx.modal_coeffs, x );

    Real P = Real( 0 );
    for( Index l = 0; l <= approx.L; ++l )
      P += approx.p_coeffs[ l ] * ReferenceElementType::legendreP( l, x );

    Real Q = Real( 0 );
    for( Index m = 0; m <= approx.M; ++m )
      Q += approx.q_coeffs[ m ] * ReferenceElementType::legendreP( m, x );

    if( std::abs( Q ) < fallback_tol_ )
      return polynomialFromModal_( approx.modal_coeffs, x );

    return P / Q;
  }

  std::vector< Real > reconstruct( const Real* nodal_vals, const std::vector< Real >& output_nodes ) const
  {
    const auto c = modalCoeffs( nodal_vals );
    const auto approx = buildApproximant( c );
    std::vector< Real > out( output_nodes.size() );
    for( std::size_t i = 0; i < output_nodes.size(); ++i )
      out[ i ] = evaluate( approx, output_nodes[ i ] );
    return out;
  }

  Field reconstruct( const Field& u ) const
  {
    Field result( u.numElements(), Np_ );
    std::vector< Real > gll_nodes( Np_ );
    for( Index i = 0; i < Np_; ++i )
      gll_nodes[ i ] = ref_.nodes()[ i ];

    for( Index k = 0; k < u.numElements(); ++k ) {
      const auto vals = reconstruct( u.elementPtr( k ), gll_nodes );
      Real* rk = result.elementPtr( k );
      for( Index i = 0; i < Np_; ++i )
        rk[ i ] = vals[ i ];
    }
    return result;
  }

private:
  EigenMat buildBMatrix_( const std::vector< Real >& c ) const
  {
    EigenMat B = EigenMat::Zero( N_ + 1, M_ + 1 );
    const Index nq = ref_.numDOF();

    for( Index q = 0; q < nq; ++q ) {
      const Real xq = ref_.nodes()[ q ];
      const Real wq = ref_.weights()[ q ];

      Real uh = Real( 0 );
      for( Index n = 0; n <= N_; ++n )
        uh += c[ n ] * ReferenceElementType::legendreP( n, xq );

      for( Index k = 0; k <= N_; ++k ) {
        const double Pk = static_cast< double >( ReferenceElementType::legendreP( k, xq ) );
        for( Index m = 0; m <= M_; ++m ) {
          const double Pm = static_cast< double >( ReferenceElementType::legendreP( m, xq ) );
          B( k, m ) += static_cast< double >( wq ) * Pk * Pm * static_cast< double >( uh );
        }
      }
    }
    return B;
  }

  Real polynomialFromModal_( const std::vector< Real >& c, Real x ) const
  {
    Real value = Real( 0 );
    for( Index n = 0; n < static_cast< Index >( c.size() ); ++n )
      value += c[ n ] * ReferenceElementType::legendreP( n, x );
    return value;
  }

  const ReferenceElementType& ref_;
  Index L_, M_, N_, Np_;
  Real fallback_tol_;
  RefMatrix vinv_;
};

} // namespace TNL::DGM
