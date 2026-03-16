#pragma once

template
<
  class RealType,
  class VectorType,
  class MatrixType,
  int number_of_elements,
  VectorType degree_of_freedom
>
class Mesh
{
public:
  // constructor
  Mesh()
  {

  }

private:
  VectorType _vertex_cord;
  MatrixType _EToV;
  VectorType _x_cord;
  int _number_elements;
  VectorType _dof;
};
