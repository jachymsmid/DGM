# DGM solver for 1D advection equation

This solver was written along the book Nodal Discontinuous Galerkin Methods by Hesthaven and Warburton.

# DGM

# Modules

## JacobiP.hpp
This script evaluates the Jacobi polynomial at a point or for all elements of a vector.
The three term recurrence formula is used for the construction.

``` math
P_{n+1}^{\alpha,\beta}(x) = [(b_n+c_n x)P_n^{\alpha,\beta}-d_n P_{n-1}^{\alpha,\beta}]/a_n \\
a_n = 2(n+1)(n+\alpha+\beta+1)(2n+\alpha+\beta), \\
b_n = (2n+\alpha+\beta+1)(\alpha^2-\beta^2), \\
c_n = (2n+\alpha+\beta)(2n+\alpha+\beta+1)(2n+\alpha+\beta+2) \\
d_n = 2(n+\alpha)(n+\beta)(2n+\alpha+\beta+2)
```
To initialize the recurrence we need to know the first two terms, they are given as
``` math
P_0^{\alpha,\beta} = 0 \\
P_1^{\alpha,\beta} = [\alpha-\beta+(\alpha+\beta+2)x]/2
```
