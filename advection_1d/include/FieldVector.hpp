/**
 * @file FieldVector.hpp
 * @brief Element-local field storage for the DG solver.
 *
 * Provides a contiguous container holding K elements each with Np degrees
 * of freedom. Supplies element pointer accessors and deep-copy utilities.
 */

#pragma once

#include <TNL/Containers/Vector.h>

namespace TNL::DGM {

/**
 * @class FieldVector
 * @brief Element-local storage for the DG solution field.
 *
 * Stores K elements each with Np degrees of freedom in a contiguous
 * TNL::Containers::Vector. Provides element pointer accessors and simple
 * utilities used throughout the solver.
 */
template
<
  class Real   = double,
  class Device = TNL::Devices::Host,
  class Index  = int
>
class FieldVector
{
public:
  using Storage = TNL::Containers::Vector<Real, Device, Index>;

  // Constructors
  /**
   * @brief Default construct an empty FieldVector.
   *
   * Results in numElements() == 0 and totalSize() == 0.
   */
  FieldVector() = default;
  /**
   * @brief Construct a FieldVector sized for K elements with Np DOFs each.
   *
   * Allocates contiguous storage and initializes values to zero.
   * @param K Number of elements
   * @param Np Degrees of freedom per element
   */
  FieldVector(Index K, Index Np) : K_(K), Np_(Np), data_(K * Np) { data_.setValue(0); }

  // Getters
  /**
   * @brief Number of mesh elements (K).
   * @return K
   */
  Index numElements() const { return K_; }
  /**
   * @brief Degrees of freedom per element (Np).
   * @return Np
   */
  Index numDOF() const { return Np_; }
  /**
   * @brief Total number of stored degrees of freedom (K * Np).
   * @return total size
   */
  Index totalSize() const { return K_ * Np_; }

  // Element-local view as pointer (host only for now)
  /**
   * @brief Host pointer to the first DOF of element k.
   * @param k element index
   * @return pointer to element-local data
   */
  Real* elementPtr(Index k) { return data_.getData() + k * Np_; }
  /**
   * @brief Const host pointer to the first DOF of element k.
   * @param k element index
   * @return const pointer to element-local data
   */
  const Real* elementPtr(Index k) const { return data_.getData() + k * Np_; }

  /**
   * @brief Access to underlying storage container.
   * @return reference to the storage vector
   */
  Storage& data() { return data_; }
  /**
   * @brief Const access to underlying storage container.
   * @return const reference to the storage vector
   */
  const Storage& data() const { return data_; }

  // Deep copy
  /**
   * @brief Deep copy from another FieldVector.
   * @param other source FieldVector to copy
   */
  void copyFrom(const FieldVector& other) { data_ = other.data_; }

private:
  Index K_{0}, Np_{0};
  Storage data_;
};

} // namespace DG
