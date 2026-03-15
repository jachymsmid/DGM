#set text(font: "New Computer Modern", lang: "en")
#show link: set text(fill: blue)
#set page(
  paper: "a4",
  number-align: center,
  numbering: "— 1 —",
)
#set heading(numbering: "1.")
#set document(
  title: [#align(center)[DG method for the 1D advection equation]],
  author: "Jachym Smid",
  date: auto
)


#title()

This text was mainly inspired by #cite(<hesthaven2008nodal>). Code supporting this text can be found on #link("https://github.com/jachymsmid/DGM/tree/main/advection_1d")[Github].

#heading(level: 1, numbering: none, "Symbols")

- $u_h$ - approximate solution obtained on a mesh with step $h$
- $(dot.op, dot.op)_(L^2)$ - scalar product in $L^2$
- $P_k (I_j)$ - space of polynomoials of up to order $k$ on an interval $I_j$
- $phi_i$ - $i$-th basis function
- $bold(f)^*$ - numerical flux
- $brace.l.double f brace.r.double = (u^- + u^+)/2$
- $bracket.l.stroked f bracket.r.stroked = hat(bold(n))^- dot.op bold(u)^- + hat(bold(n))^+ dot.op bold(u)^+$
- $hat(u)$ - vector of coefficients
- $l_i(x)$ - Lagrange polynomoial

= Problem overview
The general one dimensional advection equation is of the following form:
$ (partial bold(u))/(partial t) + (partial bold(f)(bold(u)))/(partial x) = 0 $
where $bold(f) = [f_1 (bold(u)),f_2 (bold(u)), dots, f_n (bold(u))]^T$ is the physical flux and $bold(u) = bold(u)(x,t) = [u_1, u_2, dots, u_n]^T$ is the unknown vector valued function.
This equation is hyperbolic.

For this problem to be well defined we need to impose initial and boundary conditions.
$
bold(u)(x, 0) = bold(u)_0 (x)\
cal(B)_L bold(u)(L, t) = bold(g)_1 (t) " at " x = L\
cal(B)_R bold(u)(R, t) = bold(g)_2 (t) " at " x = R
$
where the sum of the ranks of the boundary operators $cal(B)_L,cal(B)_R$ equals the number of the required inflow conditions.

== Weak formulation of the problem

To obtain the weak formulation we multiply the equation by a test function $bold(v)_h in V_h^k$ and integrate over $I^j$.
$ integral_(Omega) bold(v)^T (partial bold(u))/(partial t) dif x + integral_(Omega) bold(v)^T (partial bold(f) (bold(u)))/(partial x) dif x = 0. $
Using per-partes on the second integral we obtain
$ integral_(Omega) bold(v)^T (partial bold(u))/(partial t) dif x - integral_(Omega) (partial bold(v)^T)/(partial x) bold(f)(bold(u)) dif x + integral_(partial Omega) hat(bold(n)) dot.op bold(v)^T bold(f) (bold(u)) dif x = 0 $

= Discontinous Galerkin Method

== Spatial discretization

We split our domain (interval $I = chevron.l L, R chevron.r$) into $N$ elements $D^k= chevron.l x_(j-2/2), x_(j+1/2) chevron.r, thick j = 1,2,dots,N,$ here $x_j$ is the center of the element.
We now formulate the local weak formulation
$
integral_(D^k) bold(v)^T (partial bold(u))/(partial t) dif x - integral_(D^k) (partial bold(v)^T)/(partial x) bold(f)(bold(u)) dif x = -bracket.l bold(v)^T bold(f) (bold(u)) bracket.r_(x_(j-1/2))^(x_(j+1/2)) = integral_(partial D^k) hat(bold(n)) dot.op bold(f)^*(bold(u)) bold(v) dif x
$
We now have to define $bold(f)(bold(u)(x_(j-1/2), t))$, because the function is multivalued at the element interfaces.
To do that we define a numerical flux $bold(f)^* = bold(f)^* (bold(u)^-, bold(u)^+)$. The flux must be consistent i.e. $bold(f)^* (a,a) = bold(f) (a).$ One such numerical flux could be the local Lax-Friedrichs numerical flux.
$
bold(f)^* = brace.l.double bold(f)(bold(u)) brace.r.double = C/2 bracket.l.stroked bold(u) bracket.r.stroked
$
where the local constant $C$ is determined by the maximum eigenvalue (the spectral radius) of the flux Jacobian
$
C = max_i lambda_i = rho (bold(f)_(bold(u))) = rho ((partial bold(f))/(partial bold(u)))
$

== Basis functions

We assume that the solution $bold(u) (x,t)$ can be expressed as a direct sum of local piecewise polynomial solution of degree $K$
$
bold(u) (x,t) approx bold(u)_h (x,t) = plus.o.big_(k=1)^K bold(u)_h^k (x^k, t)
$
We define a function space
$
V_h^k = {phi in L^2(I) thick : thick phi bar.v_(D^k) in P_k (D^k), thick j = 1,2,dots,N},
$
where $P_k (D^k)$ is a space of polynomials of at most degree $k$ on the interval $I_j$.

Now we can express the local approximate solution in the basis of the space $V_h^k$
$
x in D^k quad : quad u_h^k = sum_(n)^N hat(bold(u))_n^k (t) phi_n (x) = sum_j^N bold(u)_h^k (x_i, t) l_i (x),
$
where $hat(bold(u))_n$ is a vector of coefficients and $l_i$ is the Lagrange polynomial. The first expression is said to be modal representation and the second nodal.

Following the Galerkin approach we replace the test function with each of the basis functions. This yields a system of $N+1$ equations
$
integral_(D^k) phi_i partial/(partial t) (sum_(j=0)^N hat(bold(u))_j^k phi_j) dif x - integral_(D^k) (partial phi_i)/(partial x) bold(f) (bold(u)_h^k) dif x + integral_(partial D^k) hat(bold(n)) dot.op phi_i bold(f)^* dif x = 0
$

When choosing the basis functions for the space $V_h^k$ we have many options, the easiest one is the *monomial* basis.

=== Monomial basis

For easier computation we construct our basis on a refrence element $hat(I) = chevron.l -1, 1 chevron.r$
and then we map it onto each element $I_j$ using an affine transformation. Which is given by
$ x = x_j + (1+xi)/2 h_j, $
where $xi in hat(I)$ is a local variable on $hat(I)$.
The monomial basis then is
$ hat(phi)^(m) (xi) = xi^m, quad m = 1, 2, dots k $ 
For the solver itself we need to evaluate integrals of a product of two basis functions on $I_j$, we do that by substituting onto the reference interval
$
integral_(I_j) phi_j^m phi_j^n dif x = integral_(-1)^1 hat(phi)^m hat(phi)^n h/2 dif xi =\
= h/2 integral_(-1)^1 xi^(m+n) dif xi = cases(h/(m+n+1) quad &: quad (m+n) percent 2 = 0,0 quad &: quad (m+n) percent 2 != 0)
$
These form the mass matrix $M.$

Next we will also need integrals involving a derivative of a basis function
$
integral_(I_j) phi_j^m (dif phi_j^n)/(dif x) dif x = integral_(hat(I)) hat(phi)^m (dif hat(phi)^n)/(dif xi) 2/h h/2 dif xi =\
= integral_(hat(I)) xi^m n xi^(n-1) dif xi = cases(0 quad &: quad (m+n) percent 2 = 0,(2n)/(m+n) quad &: quad (m+n)percent 2 != 0)
$
These integrals form the matrix $S$.

Next we need to know the values of the basis functions on the cell boundaries:
$
phi_j^m (x_(j-1/2)) = (-1)^m,\
phi_j^m (x_(j+1/2)) = 1
$
those values will form the matrix $B$.

The problem with this basis is that the mass matrix M is poorly conditioned. To solve this we employ the Gram-Schmidt process to orthonormalize the basis yielding us the Legender polynomial basis.

=== Legender polynomials

Legender polynomials are the result of orthonormalizing the monomial basis. They are a part of a larger family of polynomials called the Jacobi polynomials.

An easy way to calculate the polynomials is through the three term reccurence.
$
P_(n+1)^(alpha,beta)(x) &= [(b_n+c_n x)P_n^(alpha,beta)-d_n P_(n-1)^(alpha,beta)]/a_n\
&a_n = 2(n+1)(n+alpha+beta+1)(2n+alpha+beta),\
&b_n = (2n+alpha+beta+1)(alpha^2-beta^2),\
&c_n = (2n+alpha+beta)(2n+alpha+beta+1)(2n+alpha+beta+2)\
&d_n = 2(n+alpha)(n+beta)(2n+alpha+beta+2)
$
To start the reccurence we need the first two terms. They are given as
$
&P_0^(alpha,beta) = 1\
&P_1^(alpha,beta) = [alpha-beta+(alpha+beta+2)x]/2
$
Now to obtain the Legender polynomials we set the two coefficients $alpha, beta$ to zero.

This formula is described in #cite(<karniadakis2013spectral>) in the appendix A.

Using this basis the mass matrix M becomes the identity matrix (from the orthonormal character of the polynomials) and the problem of conditioning is resolved.
$
M = (phi_i, phi_j)_(L^2) = II
$

== Nodal representation

//Now we need to compute the $U_n$, which is given as
//$ U_n = (u, phi_n)_(L^2) $
//This identity follows from residual orthonormality.
//$
//(u-u_h, phi_n)_(L^2) = 0\
//(u, phi_n)_(L^2) - (u_h, phi_n)_(L^2) = 0\
//(u,phi_n)_(L^2) = (U_m phi_m, phi_n)_(L^2) = U_m (phi_m, phi_n)_(L^2) = U_m delta_(m n) = U_n\
//(u, phi_n)_(L^2) = U_n
//$
//For a general function this integral is nontrivial, we therefor approximate it with a sum. In particular we use the Gaussian quadrature.
//$
//U_n = sum_i^N u(r_i) phi_n (r_i) w_i,
//$
//where $r_i$ are suitably chosen points and $w_i$ are weights. This quadrature is exact for polynomials of order $2 N -1$.

We define $hat(bold(u))_n$ such that the approximation is interpolatory. Meaning
$
u(xi_i) = sum hat(bold(u))_n psi_n (xi_i),
$
where $xi_i$ represents distinct grid nodes. We can now write
$ bold(u) = cal(V) hat(bold(u)), $
where $cal(V)_(i j) = phi_j (xi_i)$ is the Vandermonde matrix and $u_i = u(xi_i)$. This matrix connects the modes $U$ and nodal values $u(xi_i)$. We would like the matrix to be well conditioned. We already chose the basis polynomials, so only the nodes $xi_i$ are left to define the Vandermonde matrix.

We first recognize that if
$ bold(u)(r) approx bold(u)_h (r) = sum hat(bold(u))_n phi_n $
is an interpolant ($bold(u)(xi_i) = bold(u)_h (xi_i)$), then we can write it as
$ bold(u)(r) approx bold(u)_h (r) = sum bold(u)(xi_i) l_i (r), $
where $l_i (r)$ are the Lagrange polynomials.

Without other comments the Legender-Gauss-Lobatto nodes were chosen.

== Local operators
Up till now we were developing a sensible local representation of the approximate solution. We should now discuss the various local operators.

=== Mass matrix
The local mass matrix $M^k$ is given as
$ M_(i j)^k = integral_(x_l^k)^(x_r^k) l_i^k (x) l_j^k (x) dif x = h_k/2 integral_(-1)^1 l_i (r) l_j (r) dif r = h_k/2 (l_i, l_j)_(L^2) = h_k/2 M_(i j), $
where the coefficient $h_k/2$ is a Jacobian coming from the affine transformation and $M$ is the mass matrix defined on the reference element.

We know that
$ l_i (r) = (cal(V)^T)^(-1) phi_n (r). $
From which we get
$
M_(i j) = integral_(-1)^1 (cal(V)^T)_(i n)^(-1) phi_n (r) (cal(V)^T)_(j m)^(-1) phi_m (r) dif r= (cal(V)^T)_(i m)^(-1) (cal(V)^T)_(j m)^(-1) (phi_n, phi_m)_(L^2) =\
=(cal(V)^T)_(i m)^(-1) (cal(V)^T)_(j m)^(-1)
$
Thus
$ M^k = h_k/2 M = h_k/2 (cal(V) cal(V)^T)^(-1) $

=== Stiffness matrix
The local stiffness matrix is given as
$ S_(i j)^k = integral_(x_l^k)^(x_r^k) l_i^k (x) dif/(dif x) l_j^k (x) dif x = integral_(-1)^1 l_i (r) dif/(dif r) l_j (r) dif r = S_(i j) $
To simply computate this we define the differentiation matrix.
$ D_(i j) = dif/(dif r) l_j bar.v_(r_i) $
Consider the product
$
(M D)_(i j) = M_(i n) D_(n j) = sum_n integral_(-1)^1 l_i (r) l_n (r) dif/(dif r) l_j bar.v_(r_n) = integral_(-1)^1 l_i (r) sum_n (dif l_j)/(dif r) bar.v_(r_n) l_n (r) dif r =\
= integral_(-1)^1 l_i (r) (dif l_j)/(dif r) dif r = S_(i j)\
M D = S
$
The entries of the differentiation matrix can be found directly
$ D = cal(V)_r cal(V)^(-1), $
where $cal(V)_r$ is the Vandermonde matrix assembled from differentiated Legendre polynomials. 

=== Surface integral
This operator is responsible of extracting the surface terms of the form
$ integral.cont_(-1)^1 hat(bold(n)) dot.op (u_h - u^*) l_i (r) dif r = (u_h - u^*) bar.v_(r_n) bold(e)_n - (u_h - u^*) bar.v_(r_1) bold(e)_1, $
where $bold(e)_i$ is a zero vector with $1$ at index $i$.

== Time discretization

We use explicit RK4 to integrate the semidiscrete system
$ (dif u_h)/(dif t) = cal(L)_h (u_h, t) $
$
&k_1 = cal(L)_h (u_h^n, t^n)\
&k_2 = cal(L)_h (u_h^n + 1/2 Delta t k_1, t^n + Delta t)\
&k_3 = cal(L)_h (u_h^n + 1/2 Delta t k_2, t^n + Delta t)\
&k_4 = cal(L)_h (u_h^n + 1/2 Delta t k_3, t^n + Delta t)\
&u_h^(n+1) = u_h^n + 1/6 Delta t (k_1 + k_2 + k_3 + k_4)
$

= Mesh

Unsturctured mesh was chosen as the type of spatial discretization representation.
I opted for the mesh structure from #link("https://tnl-project.org/")[TNL]


= Example problem

We are given the following problem
$
partial_u + 2 pi partial_u = 0, quad x in chevron.l 0, 2 pi chevron.r\
u(x,0) = sin(x)\
u(0,t) = -sin(2 pi t)
$

We assume the solution can be approximated as
$ u_h^k = sum_i u_h^k (x_i^k, t) l_i^k (x). $
And that the direct sum of these solutions is an approximation to the global solution. This yields the local semidiscrete scheme
$ M^k (dif u_h^k)/(dif t) + 2 pi S u_h^k = [l^k (x) (2 pi u_h^k - (2 pi u)^*)]_(x_l^k)^(x_r^k) $

#bibliography("sources.bib")
