# DGM solver for 1D advection equation

This solver was written along the book Nodal Discontinuous Galerkin Methods by Hesthaven and Warburton.

# DGM

# Modules

## JacobiP.hpp
This script evaluates the Jacobi polynomial at a point or for all elements of a vector.
The three term recurrence formula is used for the construction.

``` math
\begin{align*}
P_{n+1}^{\alpha,\beta}(x) &= [(b_n+c_n x)P_n^{\alpha,\beta}-d_n P_{n-1}^{\alpha,\beta}]/a_n \\
&a_n = 2(n+1)(n+\alpha+\beta+1)(2n+\alpha+\beta), \\
&b_n = (2n+\alpha+\beta+1)(\alpha^2-\beta^2), \\
&c_n = (2n+\alpha+\beta)(2n+\alpha+\beta+1)(2n+\alpha+\beta+2) \\
&d_n = 2(n+\alpha)(n+\beta)(2n+\alpha+\beta+2)
\end{align*}
```

To initialize the recurrence we need to know the first two terms, they are given as

``` math
\begin{align*}
&P_0^{\alpha,\beta} = 0 \\
&P_1^{\alpha,\beta} = [\alpha-\beta+(\alpha+\beta+2)x]/2
\end{align*}
```
## JacobiGassQuadrature.hpp
This script computes the nodes and weights for the Gauss quadrature.
``` math
\begin{align*}
\int_a^b f(x)\ dx \approx \sum_{i=1}^N w_i g(x_i)
\end{align*}
```
The algorithm used is the one presented in [Calculation of Gauss Quadrature Rules](https://web.stanford.edu/class/cme335/spr11/S0025-5718-69-99647-1.pdf) by Golub and Welsh. Although there have been many improvments since, both in precision and speed (see [FAST AND ACCURATE COMPUTATION OF GAUSS–LEGENDRE
AND GAUSS–JACOBI QUADRATURE NODES AND WEIGHTS](https://epubs.siam.org/doi/10.1137/120889873)), I opted for the simpler GW algorithm.


