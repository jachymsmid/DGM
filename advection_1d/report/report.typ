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

#heading(level: 1, numbering: none, "Symbols and definitions")

- $u_h$ - approximate solution obtained on a mesh with step $h$
- $(u, v)_(L^2(Omega)) = integral_Omega u v dif x$ - scalar product in $L^2$
- $bar.v.double u bar.v.double _(L^2 (Omega))^2 = (u, u)_(L^2 (Omega))$
- $P_j (D^k)$ - space of polynomoials of up to order $j$ on an interval $D^k$
- $phi_i, psi_i$ - $i$-th basis function
- $u^-$ - limit of $u$ when approaching the boundary of an element from the left
- $u^+$ - limit of $u$ when approaching the boundary of an element from the right
- $f^* = f^* (u^-, u^+)$ - numerical flux
- $hat(bold(n))$ - outward unit normal
- $brace.l.stroked u brace.r.stroked = (u^- + u^+)/2$
- $bracket.l.stroked u bracket.r.stroked = hat(bold(n))^- u^- + hat(bold(n))^+ u^+$
- $bracket.l.stroked bold(u) bracket.r.stroked = hat(bold(n))^- dot.op bold(u)^- + hat(bold(n))^+ dot.op bold(u)^+$
- $l_i (x)$ - $i$-th Lagrange polynomoial
- $N$ - order of polynomial approximation
- $N_p = N + 1$ - number of points in each element (number of degrees of freedom)

= Problem overview
The general one dimensional advection equation is of the following form:
$ (partial bold(u))/(partial t) + (partial bold(f)(bold(u)))/(partial x) = bold(s)(x,t) quad [x,t] in chevron.l L, R chevron.r times RR^+ $
where $bold(f) = [f_1 (bold(u)),f_2 (bold(u)), dots, f_n (bold(u))]^T$ is the physical flux, $bold(u) = bold(u)(x,t) = [u_1, u_2, dots, u_n]^T$ is the unknown vector valued function and $bold(s)(x,t)$ is the source term.
This equation is hyperbolic.

For this problem to be well defined we need to impose initial and boundary conditions.
$
bold(u)(x, 0) = bold(u)_0 (x)\
cal(B)_L bold(u)(L, t) = bold(g)_1 (t) " at " x = L\
cal(B)_R bold(u)(R, t) = bold(g)_2 (t) " at " x = R
$
where the sum of the ranks of the boundary operators $cal(B)_L,cal(B)_R$ equals the number of the required inflow conditions.

== Weak formulation of the problem

To obtain the weak formulation we multiply the equation by a test function $bold(v) in V$ and integrate over the domain $Omega$.
$ integral_(Omega) bold(v)^T (partial bold(u))/(partial t) dif x + integral_(Omega) bold(v)^T (partial bold(f) (bold(u)))/(partial x) dif x = integral_(Omega) bold(v)^T bold(s)(x,t) dif x. $
Using per-partes on the second integral we obtain
$ integral_(Omega) bold(v)^T (partial bold(u))/(partial t) dif x - integral_(Omega) (partial bold(v)^T)/(partial x) bold(f)(bold(u)) dif x + integral_(partial Omega) hat(bold(n)) dot.op bold(v)^T bold(f) (bold(u)) dif x = integral_(Omega) bold(v)^T bold(s)(x,t) dif x, $
$hat(bold(n))$ here is a unit outward normal. Notice that we express the term $[bold(v)^T bold(f)(bold(u))]_L^R$ in a integral form $integral_(partial Omega) hat(bold(n)) dot.op bold(v)^T bold(f) (bold(u)) dif x$, this will later help us when generalazing the scheme.

= Discontinous Galerkin Method

== Spatial discretization

We split our domain ($Omega = chevron.l L, R chevron.r$) into $N$ elements $D^j= chevron.l x_(j-1/2), x_(j+1/2) chevron.r, thick j = 1,2,dots,N,$ here $x_j$ is the center of the element and $x_(1/2) = L$, $x_(N+1/2) = R$.
$ Omega approx Omega_h = limits(union.big)_(k=1)^K D^k $
We now formulate the local weak formulation
$
integral_(D^k) bold(v)^T (partial bold(u))/(partial t) dif x = integral_(D^k) (partial bold(v)^T)/(partial x) bold(f)(bold(u)) dif x - integral_(partial D^k) hat(bold(n)) dot.op bold(v)^T bold(f)(bold(u)) dif x + integral_(D^k) bold(v)^T bold(s)(x,t) dif x
$
But now we have a problem, because $u$ is double valued at the boundaries. To solve this we define a numerical flux $bold(f)^* = bold(f)^* (bold(u)^-, bold(u)^+)$. The equation then becomes
$
integral_(D^k) bold(v)^T (partial bold(u))/(partial t) dif x = integral_(D^k) (partial bold(v)^T)/(partial x) bold(f)(bold(u)) dif x - integral_(partial D^k) hat(bold(n)) dot.op bold(v)^T bold(f)^*(bold(u)) dif x + integral_(D^k) bold(v)^T bold(s)(x,t) dif x
$
The flux must be consistent i.e. $bold(f)^* (a,a) = bold(f) (a).$ One such numerical flux could be the local Lax-Friedrichs numerical flux.
$
bold(f)^* = brace.l.stroked bold(f)(bold(u)) brace.r.stroked + C/2 bracket.l.stroked bold(u) bracket.r.stroked
$
where the local constant $C$ is determined by the maximum eigenvalue (the spectral radius) of the physical flux Jacobi matrix.
$
C = max_i lambda_i = rho (bold(f)_(bold(u))) = rho ((partial bold(f))/(partial bold(u)))
$

== Basis functions

We assume that the solution $bold(u) (x,t)$ can be expressed as a direct sum of local piecewise polynomial solutions
$
bold(u) (x,t) approx bold(u)_h (x,t) = plus.o.big_(k=1)^K bold(u)_h^k (x^k, t)
$
We define a local function space $V_h^k$ such that
$
V_h = plus.o.big_(k=1)^K V_h^k\
V_h^k = {phi in L^2(D^k) thick : thick phi bar.v_(D^k) in P_j (D^k), thick j = 1,2,dots,N},
$
where $P_j (D^k)$ is a space of polynomials of at most degree $k$ on the interval $D^j.$

Now we can express the local approximate solution in the basis of the space $V_h^k$
$
x in D^k quad : quad bold(u)_h^k = sum_(n)^(N_p) hat(bold(u))_n^k (t) phi_n (x) = sum_i^(N_p) bold(u)_h^k (x_i, t) l_i (x),
$
where $hat(bold(u))_n$ is a vector of coefficients and $l_i$ is the $i$-th interpolating Lagrange polynomial and $x_i$ are distinct. The first expression is said to be modal representation and the second nodal. We won't discuss the modal formulation in this text.

#line(length: 100%)
*NOTE:* only linear problems without source term for now

*TODO:* theory for nonlinear problems, add source term
#line(length: 100%)

Following the Galerkin approach we replace the test function with each of the basis functions. This yields a system of $N+1$ equations for each element
$
integral_(D^k) l_i (x) (partial u_h (x_j, t) l_j (x))/(partial t) dif x = integral_(D^k) a u_h (x_i, t) l_i (x) (dif l_j (x))/(dif x) dif x - integral_(partial D^k) l_i f^* dif x\
sum_(j=1)^N_p u_h (x_j, t) integral_(D^k) l_i (x) (partial l_j (x))/(partial t) dif x = a integral_(D^k) l_i (dif l_j (x))/(dif x) dif x - integral_(partial D^k) l_i f^* dif x\
i = 1,2,dots,N+1\
$

When choosing the basis of the space $V_h^k$ we have many options, but only one makes sense, the Legendre polynomials, $psi_n$ is a Legendre polynomial of order $n-1$. It will come in handy that they are orthogonal, and we will also normalize them in the $L^2$ norm, but more on them later.

=== Legender polynomials
... blah blah blah ...

An easy way to calculate the polynomials is through the three term reccurence using the Bonnet's formula [cite]
$
P_n = ((2 n + 1) x P_(n-1) (x) - n P_(n-2) (x)) / (n + 1);\
$
To start the reccurence we need the first two terms. They are given as
$
P_0 (x) = 1 quad P_1 (x) = x
$
Now to obtain the normalized (in the $L^2$ norm) Legender polynomials $tilde(P)_n (x)$ we multiply each polynomial by an appropriate coefficient
$
tilde(P)_n (x) = sqrt((2 n + 1)/2) P_n (x)
$

To obtain a derivative of the Legendre polynomial we use the formula [cite]
$
P_n^' (x) = (n (P_(n-1) (x) - x P_(n) (x)))/(1-x^2)
$
Notice that the formula doesn't hold for $x = plus.minus 1$, for that we use a the following property
$
P_n^' (1) =  (n (n+1))/2
$
and
$
P_n^' (-x) = (-1)^(2 n - 1) P_n^' (x)
$

#heading(level: 4, numbering: none, "Important properties")
+ $ P_n (-x) = (-1)^n P_n (x) $
+ $ integral_(-1)^1 P_n (x) P_m (x) dif x = 0 " for " n != m $
+ $ integral_(-1)^1 P_n (x)  dif x = 0 " for " n gt.eq 1 $
+ ...

== Nodal representation

Let's now discuss the nodal representation in greater detail. We define $hat(bold(u))_n$ such that the approximation is interpolatory i.e.
$
u(x_i) = sum hat(u)_n psi_n (x_i),
$
where $xi_i$ represents distinct grid nodes. We can now write
$ bold(u) = cal(V) hat(bold(u)), $
where $cal(V)_(i j) = phi_j (xi_i)$ is the Vandermonde matrix and $u_i = u(xi_i)$. This matrix connects the modes $U$ and nodal values $u(xi_i)$. We would like the matrix to be well conditioned. We already chose the basis polynomials, so only the nodes $xi_i$ are left to define the Vandermonde matrix.

We first recognize that if
$ bold(u)(r) approx bold(u)_h (r) = sum hat(bold(u))_n phi_n $
is an interpolant ($bold(u)(xi_i) = bold(u)_h (xi_i)$), then we can write it as
$ bold(u)(r) approx bold(u)_h (r) = sum bold(u)(xi_i) l_i (r), $
where $l_i (r)$ are the Lagrange polynomials defined as
$
l_i (x) = limits(Pi)_(j=1, j != i)^N_p thick (x-x_i)/(x_i - x_j)
$
An important property of the Lagrange polynomial is that
$
l_i (x_j) = delta_(i j),
$
where $delta_(i j)$ is the Kronecker delta.

Without other comments the Legender-Gauss-Lobatto nodes were chosen.

== Local operators
Up till now we were developing a sensible local representation of the approximate solution. We should now discuss the various local operators in the nodal formulation of DGM.

=== Mass matrix
The local mass matrix $M^k$ is given as
$ M_(i j)^k = integral_(x_l^k)^(x_r^k) l_i^k (x) l_j^k (x) dif x = h_k/2 integral_(-1)^1 l_i (r) l_j (r) dif r = h_k/2 (l_i, l_j)_(L^2) = h_k/2 M_(i j), $
where the coefficient $h_k/2$ is a Jacobian coming from the affine transformation and $M$ is the mass matrix defined on the reference element.

We know that
$ l_i (r) = (cal(V)^T)^(-1)_(i n) phi_n (r). $
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
$ integral_(partial D^k) hat(bold(n)) f^* l_i (r) dif x =  f^* bar.v_(x_R) bold(e)_n - u^* bar.v_(x_L) bold(e)_1, $
where $bold(e)_i$ is a zero vector with $1$ at index $i$. We can rewrite this using a $N_p times 2$ matrix $cal(E)$ that is zero except for upper right corener that is $1$ and the lower left corner that is $-1$ like so
$ integral_(partial D^k) hat(bold(n)) f^* l_i (r) dif x = cal(E) dot.op (f^*_R, f^*_L)^T $

== Time discretization

Two explicit RK4 mmethods were implemented to integrate the semidiscrete system
$ (dif u_h)/(dif t) = cal(L)_h (u_h, t) $

=== ERK
General explicit four step Runge-Kutta method.
$
&k_1 = cal(L)_h (u_h^n, t^n)\
&k_2 = cal(L)_h (u_h^n + 1/2 Delta t k_1, t^n + Delta t)\
&k_3 = cal(L)_h (u_h^n + 1/2 Delta t k_2, t^n + Delta t)\
&k_4 = cal(L)_h (u_h^n + 1/2 Delta t k_3, t^n + Delta t)\
&u_h^(n+1) = u_h^n + 1/6 Delta t (k_1 + k_2 + k_3 + k_4)
$

=== LSERK
Low storage explicit five step Runge-Kutta method.

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
$
M^k (dif u_h^k)/(dif t) - 2 pi S^T u_h^k = - [l^k (x)  (2 pi u)^*)]_(x_l^k)^(x_r^k) = - cal(E) ((2 pi u)^*, (2 pi u)^*)^T\
J^k M (dif u_h^k)/(dif t) - 2 pi S^T u_h^k = - cal(E) ((2 pi u)^*, (2 pi u)^*)^T\
(dif u_h^k)/(dif t) = 1/J^k (2 pi M^(-1) S^T u_h^k - M^(-1) cal(E) (f^*_R, f^*_L)^T)
$

We can employ some tricks from linear algebra so we don't have to invert the mass matrix, but I will not be describing them here.

... some figures ...

#bibliography("sources.bib")
