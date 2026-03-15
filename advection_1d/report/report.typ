#set text(font: "New Computer Modern", lang: "en")
#set heading(numbering: "1.")
#set document(
  title: [#align(center)[GD method for the 1D linear advection equation]],
  author: "Jachym Smid",
  date: auto
)


#title()

This text is mainly inspired by #cite(<hesthaven2008nodal>).

= TODO:
- [ ] more general formulation
- [ ] unstructured mesh
- [ ] nonlinear equation
- [ ] system of equations
- [ ] nonlinear systems

= Symbols

- $u_h$ - approximate solution obtained on a mesh with step $h$
- $(dot.op, dot.op)_(L^2)$ - scalar product in $L^2$
- $P_k (I_j)$ - space of polynomoials of up to order $k$ on an interval $I_j$
- $phi_i$ - $i$-th basis function

= Problem overview
The linear advection equation is of the following form:
$ (partial u)/(partial t) + a (partial u)/(partial x) = 0 $
where $a$ is the advection speed and $u = u(x,t)$ is the unknown scalar function.
This equation is hyperbolic.

For the initial conditions problem we have to specify the initial conditions $u(x,t) bar.v_(t = 0) = u_0 (x)$, this can be solved analytically using the method of charcteristics.
We will be solving this problem on an interval $I = chevron.l 0, L chevron.r$, with the boundary condition on the left border (assuming $a>0$) $u(x,t) bar.v_(x = 0) = alpha(t)$.

= Discontinous Galerkin Method

== Spatial discretization

We split our domain (interval $I$) into $N$ elements $I_j = chevron.l x_(j-2/2), x_(j+1/2) chevron.r, thick j = 1,2,dots,N,$ here $x_j$ is the center of the element. We will use constant spatial step $h = L/N$,
that means that $x_i = i h$. In each element we will be trying to approximate the solution with a polynomial of degree $k.$

We define a function space
$
V_h^k = {v in L^2(I) thick : thick v bar.v_(I_j) in P_k (I_j), thick j = 1,2,dots,N},
$
where $P_k (I_j)$ is a space of polynomials of at most degree $k$ on the interval $I_j$.

== Weak formulation of the problem

To obtain the weak formulation we multiply the equation by a test function $v in V_h^k$ and integrate over the interval $I_j$.
$ integral_(I_j) (partial u_h)/(partial t) v dif x + integral_(I_j) (partial u_h)/(partial x) v dif x = 0. $
Using per-partes on the second integral we obtain
$ integral_(I_j) (partial u_h)/(partial t) v dif x - integral_(I_j) (partial v)/(partial x) u_h dif x + bracket.l a u_h v bracket.r_(x_(j - 1/2))^(x_(j + 1/2)) = 0 $
Now we have to define $u_h (x_(i plus 1/2))$, to do that we define suitable numerical flux.
$ u_h (x_(x plus 1/2)) = hat(u)_h (x_(x plus 1/2)) $
One such numerical flux could be the well known upwind ($a > 0$)
$ hat(u) (x_(x plus 1/2)) = u_h^- (x_(x plus 1/2)) = limits(lim)_(x arrow x_(j + 1/2)) u_h(x) $

The integral equation then becomes
$ integral_(I_j) (partial u_h)/(partial t) v dif x - integral_(I_j) (partial v)/(partial x) u_h dif x + [ a hat(u)_h (x_(j - 1/2)) v (x_(j - 1/2)) - a hat(u)_h (x_(j + 1/2)) v (x_(j + 1/2))] = 0 $
If we sum over all the elements we get the following equation
$ sum_(j = 1)^N (integral_(I_j) (partial u_h)/(partial t) v dif x - integral_(I_j) (partial v)/(partial x) u_h dif x + [ a hat(u)_h (x_(j - 1/2)) v (x_(j - 1/2)) - a hat(u)_h (x_(j + 1/2)) v (x_(j + 1/2))]) = 0 $

== Basis functions

We assume the local solutions to be of the form
$
x in I_j quad : quad u_h^j = sum_(n)^N U_n^j (t) phi_n (x) = sum_j^N u_h^k (x_i, t) l_i (x),
$
where the first expression is said to be modal representation and the second nodal.

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

Now we need to compute the $U_n$, which is given as
$ U_n = (u, phi_n)_(L^2) $
This identity follows from residual orthonormality.
$
(u-u_h, phi_n)_(L^2) = 0\
(u, phi_n)_(L^2) - (u_h, phi_n)_(L^2) = 0\
(u,phi_n)_(L^2) = (U_m phi_m, phi_n)_(L^2) = U_m (phi_m, phi_n)_(L^2) = U_m delta_(m n) = U_n\
(u, phi_n)_(L^2) = U_n
$
For a general function this integral is nontrivial, we therefor approximate it with a sum. In particular we use the Gaussian quadrature.
$
U_n = sum_i^N u(r_i) phi_n (r_i) w_i,
$
where $r_i$ are suitably chosen points and $w_i$ are weights. This quadrature is exact for polynomials of order $2 N -1$.

But for later generalization we define $U_n$ such that the approximation is interpolatory. Meaning
$
u(xi_i) = sum U_n psi_n (xi_i),
$
where $xi_i$ represents distinct grid nodes. We can now write
$ u = cal(V) U, $
where $cal(V)_(i j) = phi_j (xi_i)$ is the Vandermonde matrix and $u_i = u(xi_i)$. This matrix connects the modes $U$ and nodal values $u(xi_i)$. We would like the matrix to be well conditioned. We already chose the basis polynomials, so only the nodes $xi_i$ are left to define the Vandermonde matrix.

We first recognize that if
$ u(r) approx u_h(r) j= sum U_n phi_n $
is an interpolant ($u(xi_i)=u_h (xi_i)$), then we can write it as
$ u(r) approx u_h(r) = sum u(xi_i) l_i (r), $
where $l_i(r)$ are the Lagrange polynomials.

Without other comments the Legender-Gauss-Lobatto nodes will be chosen.

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

Unsturctured mesh was chosen as the type of spatial discretization.


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
