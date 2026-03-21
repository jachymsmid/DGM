#pragma once

#include <TNL/Containers/Vector.h>

namespace DG {

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
  FieldVector() = default;
  FieldVector(Index K, Index Np) : K_(K), Np_(Np), data_(K * Np) { data_.setValue(0); }

  // Getters
  Index numElements() const { return K_; }
  Index numDOF() const { return Np_; }
  Index totalSize() const { return K_ * Np_; }

  // Element-local view as pointer (host only for now)
  Real* elementPtr(Index k) { return data_.getData() + k * Np_; }
  const Real* elementPtr(Index k) const { return data_.getData() + k * Np_; }

  Storage& data() { return data_; }
  const Storage& data() const { return data_; }

  // Deep copy
  void copyFrom(const FieldVector& other) { data_ = other.data_; }

private:
  Index K_{0}, Np_{0};
  Storage data_;
};

} // namespace DG
