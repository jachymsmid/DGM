#set text(font: "New Computer Modern", lang: "en")
#show link: set text(fill: blue)
#set page(
  paper: "a4",
  number-align: center,
  numbering: "1",
  footer: context {
  let num = page.numbering
  if num != none {
    set align(page.number-align.x)
    "— " + numbering(num, ..counter(page).get()) + " —"
  }}

)
#set heading(numbering: "1.")
#set document(
  title: [#align(center)[Nodal DG method for the 1D advection equation]],
  author: "Jachym Smid",
  date: auto
)

#title()

// add some interesting picture

#pagebreak()

This text was mainly inspired by #cite(<hesthaven2008nodal>). My own code supporting this text can be found on #link("https://github.com/jachymsmid/DGM/tree/main/advection_1d")[Github].

#outline()

#pagebreak()

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

#pagebreak()

= Problem overview

The general one dimensional advection equation is of the following form:
$
(partial bold(u))/(partial t) + (partial bold(f)(bold(u)))/(partial x) =
bold(s)(x,t) quad [x,t] in chevron.l L, R chevron.r times RR^+
$
where $bold(f) = [f_1 (bold(u)),f_2 (bold(u)), dots, f_n (bold(u))]^T$
is the physical flux, $bold(u) = bold(u)(x,t) = [u_1, u_2, dots, u_n]^T$
is the unknown vector valued function and $bold(s)(x,t) = [s_1, s_2, dots, s_n]$
is the source term. This equation is hyperbolic.

For this problem to be well defined we need to impose initial and boundary conditions.
$
bold(u)(x, 0) = bold(u)_0 (x)\
cal(B)_L bold(u)(L, t) = bold(g)_1 (t) " at " x = L\
cal(B)_R bold(u)(R, t) = bold(g)_2 (t) " at " x = R
$
where the sum of the ranks of the boundary operators $cal(B)_L,cal(B)_R$
equals the number of the required inflow conditions.

== Weak formulation of the problem

To obtain the weak formulation we multiply the equation by a test function $bold(v) in V$ and integrate over the domain $Omega$.
$
integral_(Omega) bold(v)^T (partial bold(u))/(partial t) dif x + integral_(Omega) bold(v)^T (partial bold(f) (bold(u)))/(partial x) dif x = integral_(Omega) bold(v)^T bold(s)(x,t) dif x.
$
Using per-partes on the second integral we obtain
$
integral_(Omega) bold(v)^T (partial bold(u))/(partial t) dif x - integral_(Omega) (partial bold(v)^T)/(partial x) bold(f)(bold(u)) dif x + integral_(partial Omega) hat(bold(n)) dot.op bold(v)^T bold(f) (bold(u)) dif x = integral_(Omega) bold(v)^T bold(s)(x,t) dif x,
$
where $hat(bold(n))$ is a unit outward normal. Notice that we express the term
$[bold(v)^T bold(f)(bold(u))]_L^R$ in a integral form
$integral_(partial Omega) hat(bold(n)) dot.op bold(v)^T bold(f) (bold(u)) dif x$,
this will later help us when moving to higher dimensions.

= Discontinous Galerkin Method

Let's focus on linear scalar advection equation without source term for now and
generalize the scheme later.

== Spatial discretization

We split our domain $Omega = chevron.l L, R chevron.r$ into $N$ elements
$
D^k= chevron.l x_(k-1/2), x_(k+1/2) chevron.r, thick k = 1,2,dots,N,
$
here $x_k$ is the center of the  $k$-th element and $x_(1/2) = L$, $x_(N+1/2) = R$.
$
Omega approx Omega_h = limits(union.big)_(k=1)^K D^k
$
We now formulate the local weak formulation
$
integral_(D^k) v (partial u)/(partial t) dif x =
integral_(D^k) (partial v)/(partial x) f(u) dif x -
integral_(partial D^k) hat(bold(n)) v f(u) dif x
$

But now we encounter a problem, because $u$ is double valued at the element boundaries.
What vallue should we choose when evaluating the $integral_(partial D^k) hat(bold(n)) v f(u) dif x$ term?
To solve this we define a numerical flux $f^* = f^* (u^-, u^+)$ that takes in both values and returns only one.
The equation then becomes
$
integral_(D^k) v (partial u)/(partial t) dif x = integral_(D^k) (partial v)/(partial x) f(u) dif x - integral_(partial D^k) hat(bold(n)) v f^*(u) dif x // + integral_(D^k) v (s)(x,t) dif x
$
The flux must be consistent i.e. $f^* (a,a) = f (a).$

#figure(
  image("img/elements_multivalue.pdf", width: 60%),
  caption: [ $u_h$ is multivalued at the elemental interfaces ],
)

*Lax-Friedrichs flux*

One such numerical flux is the local Lax-Friedrichs numerical flux.
$
f^* (u^-, u^+) = (f(u^-) + f(u^+))/2 + C/2 hat(bold(n)) (u^- - u^+)
$
where the local constant $C$ is determined by the maximum local advection speed //eigenvalue (the spectral radius) of the physical flux Jacobi matrix.
$
C = max_(u^- med lt.eq med s med lt.eq med u^+) f_u (s) //= rho (f_(u) (s)) = rho ((partial f (s))/(partial u))
$
and $hat(bold(n))$ is local unit outward normal. In the 1D example it becomes $-1$
on the left face of an element and $1$ on the right face.

*Upwind flux*

$
f^* (u^-, u^+) = cases(f(u^-) quad &: quad hat(bold(n)) C > 0,f(u^+) quad &: quad hat(bold(n)) C < 0)
$
where $C = f(u^-)$ i.e. the local advection speed.

*Godunov flux*
$
f^* (u^-, u^+) = cases(
  min(f(u^-), f(u^+)) quad &: quad hat(bold(n)) u^- < hat(bold(n)) u^+,
  max(f(u^-), f(u^+)) quad &: quad hat(bold(n)) u^- gt.eq hat(bold(n)) u^+
)
$

// *Roe flux*

== Basis functions

To approximate the solution we assume that it can be expressed as a
direct sum of local piecewise polynomial solutions, local here refers to elements
$
u (x,t) approx u_h (x,t) = plus.o.big_(k=1)^K u_h^k (x^k, t).
$
We define a local function space $V_h^k$ such that
$
V_h = plus.o.big_(k=1)^K V_h^k\
V_h^k = {
  phi in L^2(D^k) thick : thick phi bar.v_(D^k) in P_j (D^k), thick j = 1,2,dots,N
        },
$
where $P_j (D^k)$ is a space of polynomials of at most degree $j$ on the interval $D^k.$

Now we can express the local approximate solution in the basis of the space $V_h^k$
$
x in D^k quad : quad u_h^k = sum_(n)^(N_p) hat(u)_n^k (t) phi_n (x) =
sum_i^(N_p) u_h^k (x_i, t) l_i (x),
$
where $hat(u)_n$ is a vector of coefficients and $l_i$ is the $i$-th
interpolating Lagrange polynomial, we assume that $x_i in D^k$ are distinct.
The first expression is said to be modal representation and the second nodal.
We won't discuss the modal formulation any further in this text,
but it will come in handy when implementing limiters later when formulating the
scheme for general non-linear problems.

Following the Galerkin approach we replace the test function with each of
the basis functions. This yields a system of $N+1$ equations for each element
$
integral_(D^k) l_i (x) partial/(partial t)(sum_(j=1)^N_p u_h (x_j, t)) dif x = integral_(D^k) (dif l_i (x))/(dif x) sum_(j=1)^N_p a u_h (x_j, t) l_j (x) dif x - integral_(partial D^k) l_i f^* dif x\
dif/(dif t) sum_(j=1)^N_p u_h (x_j, t) integral_(D^k) l_i (x) l_j (x) dif x = sum_(j=1)^N_p a u_h (x_j, t) integral_(D^k) (dif l_i (x))/(dif x) l_j (x) dif x - integral_(partial D^k) l_i f^* dif x\
i = 1,2,dots,N+1\
$
We can rewrite this equation into it's matrix form
$
M^k (dif bold(u)_h^k)/(dif t) = a (S^k)^T bold(u)_h^k - cal(E) bold(f)^k_*,
$
where
$
&M_(i j)^k = integral_(D^k) l_i (x) l_j (x) dif x " - the mass matrix"\
&S_(i j)^k = integral_(D^k) l_i (x) (dif l_j (x))/(dif x) dif x " - the stifness matrix"\
&cal(E)_(i j) = cases(&1 quad : quad "upper left corner and lower right corner",
                &0 quad : quad "otherwise")\
&(bold(u)_h)_j^k = u_h (x_j^k, t) " - the vector of unknowns"\
&(bold(f)_h)_j^k " - the physical flux vector"\
&bold(f)^k_* = (bold(n)_L f^*(x_L^k), bold(n)_R f^*(x_R^k))^T "- numerical flux at endpoints"\
&bold(n)_L, bold(n)_R "- unit normal at the left, respectively right, element boundary"
$
This matrix equation is the semi-discrete form of the PDE.
We will discuss all the local operatros in more detail later.

When choosing the basis of the space $V_h^k$ we have many options.
But for the mass matrix to be well conditioned we will choose the Legendre polynomials.

=== Legender polynomials
Legendre polynomials are a complete set of orthogonal polynomials defined on
$chevron.l -1, 1 chevron.r$

An easy way to calculate the polynomials is through the three term reccurence
using the Bonnet's formula [cite]
$
P_n = ((2 n + 1) x P_(n-1) (x) - n P_(n-2) (x)) / (n + 1);\
$
To start the reccurence we need the first two terms. They are given as
$
P_0 (x) = 1, quad P_1 (x) = x
$
Now to obtain the normalized (in the $L^2$ norm) Legender polynomials
$tilde(P)_n (x)$ we multiply each polynomial by an appropriate coefficient
$
tilde(P)_n (x) = sqrt((2 n + 1)/2) P_n (x)
$

To obtain a derivative of the Legendre polynomial we use the formula [cite]
$
P_n^' (x) = (n (P_(n-1) (x) - x P_(n) (x)))/(1-x^2)
$
Notice that the formula doesn't hold for $x = plus.minus 1$, for that we use a the property
$
P_n^' (1) =  (n (n+1))/2
$
and
$
P_n^' (-x) = (-1)^(2 n - 1) P_n^' (x)
$

// #heading(level: 4, numbering: none, "Important properties")
// + $ P_n (-x) = (-1)^n P_n (x) $
// + $ integral_(-1)^1 P_n (x) P_m (x) dif x = 0 " for " n != m $
// + $ integral_(-1)^1 P_n (x)  dif x = 0 " for " n gt.eq 1 $

== Nodal representation

Let's now discuss the nodal representation in greater detail.
We define $hat(bold(u))_n$ such that the approximation is interpolatory i.e.
$
u(x_i) = sum hat(u)_n psi_n (x_i),
$
where $x_i$ represents distinct grid nodes. We can now write
$ bold(u) = cal(V) hat(bold(u)), $
where $cal(V)_(i j) = psi_j (x_i)$ is the Vandermonde matrix and
$u_i = u(x_i)$. This matrix connects the modes $hat(u)$ and nodal values $u(x_i)$.
We would like the matrix to be well conditioned.
We already chose the basis polynomials, so only the nodes $x_i$ are left,
to define the Vandermonde matrix.
The node's that acomplish this are the solution of
$
(1-x^2) dif/(dif x) P_N (x) = 0
$
these are known as the Legendre-Gauss-Lobatto nodes.
The LGL nodes will be noted by the greek letter $xi$.
There are $N+1$ LGL nodes for solution approximation of order $N$.

// Now we recognize that if
// $ bold(u)(r) approx bold(u)_h (r) = sum hat(bold(u))_n phi_n $
// is an interpolant ($bold(u)(xi_i) = bold(u)_h (xi_i)$), then we can write it as
// $ bold(u)(r) approx bold(u)_h (r) = sum bold(u)(xi_i) l_i (r), $
// where $l_i (r)$ are the Lagrange polynomials defined as
// $
// l_i (x) = limits(Pi)_(j=1, j != i)^N_p thick (x-x_i)/(x_i - x_j)
// $
// An important property of the Lagrange polynomial is that
// $
// l_i (x_j) = delta_(i j),
// $
// where $delta_(i j)$ is the Kronecker delta.
//
// Without other comments the Legender-Gauss-Lobatto nodes were chosen.

== Local operators
Up till now we were developing a sensible local representation
of the approximate solution. We should now discuss the various
local operators in the nodal formulation of DGM.

=== Mass matrix
The local mass matrix $M^k$ (on $k$-th element) is given as
$
M_(i j)^k = integral_(x_l^k)^(x_r^k) l_i^k (x) l_j^k (x) dif x = h_k/2 integral_(-1)^1 l_i (r) l_j (r) dif r = h_k/2 M_(i j),
$
where the coefficient $h_k/2$ is a Jacobian coming from the affine transformation and $M$ is the mass matrix defined on the reference element.

We know that
$
l_i (r) = cal(W)_(i n) phi_n (r),
$
where $cal(W) = (cal(V)^T)^(-1)$
From which we get
$
M_(i j) = integral_(-1)^1 cal(W)_(i n) phi_n (r) cal(W)_(j m) phi_m (r) dif r =
cal(W)_(i m) cal(W)_(j m) (phi_n, phi_m)_(L^2) =\
= cal(W)_(i m) cal(W)_(j m)
$
Thus
$
M^k = h_k/2 M = h_k/2 (cal(V) cal(V)^T)^(-1)
$

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
$
(cal(V)_r)_(i j) = psi'_j (xi_i)
$

But in the weak formulation $S^T$ is involved, for that we use the identity
$
S^T = D^T M
$
To recover the semi-discrete formulation we need to multiply this by the mass matrix from the left
$
M^(-1) S^T = M^(-1) D^T M = cal(V) cal(V)_r^T M = D_w
$
We call this matrix $D_w$, differentiation matrix for the weak form.

=== Surface integral
This operator is responsible of extracting the surface terms of the form
$
integral_(partial D^k) hat(bold(n)) f^* l_i (r) dif x =  (hat(bold(n)) f^*) bar.v_(x_R) bold(e)_n + (hat(bold(n)) f^*) bar.v_(x_L) bold(e)_1,
$
where $bold(e)_i$ is a zero vector with $1$ at index $i$. We can rewrite this using a $N_p times 2$ matrix $cal(E)$ that is zero except for upper right corener that is $1$ and the lower left corner that is $-1$ like so
$
integral_(partial D^k) hat(bold(n)) f^* l_i (r) dif x = cal(E) dot.op (f^*_R, f^*_L)^T
$
Again, to recover the semi-discrete form we multiply by the inverse mass matrix from the left. We call this matrix the LIFT.
$
"LIFT" = M^(-1) cal(E) = cal(V) cal(V)^T cal(E)
$

== Time discretization

Two explicit RK4 methods were implemented to integrate the semidiscrete system
$
(dif u_h)/(dif t) = cal(L)_h (u_h, t) // = 1/J^k (D_w bold(f)^k - "LIFT" bold(f)^k_*)
$

=== ERK
General explicit four stage Runge-Kutta method.
$
&k_1 = cal(L)_h (u_h^n, t^n)\
&k_2 = cal(L)_h (u_h^n + 1/2 Delta t k_1, t^n + Delta t)\
&k_3 = cal(L)_h (u_h^n + 1/2 Delta t k_2, t^n + Delta t)\
&k_4 = cal(L)_h (u_h^n + 1/2 Delta t k_3, t^n + Delta t)\
&u_h^(n+1) = u_h^n + 1/6 Delta t (k_1 + k_2 + k_3 + k_4)
$

=== LSERK
Low storage explicit five stage fourth order Runge-Kutta method. [cite the nasa paper]

=== SSPRK
strong stability preserving Runge-Kutta fourth order method [cite]
$
&k_1 = u^n + 1/2 Delta t cal(L)_h (u^n, t_n)\
&k_2 = k_1 + 1/2 Delta t cal(L)_h (k_1, t_n + 1/2 Delta t)\
&k_3 = 2/3 u^n + 1/2 k_2 + 1/6 cal(L)_h (k_2, t^n + Delta t)\
&u_h^(n+1) = k_3 + 1/2 Delta t cal(L)_h (k_3, t^n + 1/2 Delta t)
$

== Example linear problem

We are given the following problem
$
partial_t u + partial_x u = 0, quad x in chevron.l -1, 1 chevron.r\
$
with periodic boundary conditions.

We assume the solution can be approximated as
$
u_h^k = sum_i u_h^k (x_i^k, t) l_i^k (x).
$
And that the direct sum of these solutions is an approximation to the global solution. This yields the local semidiscrete scheme
$
M^k (dif u_h^k)/(dif t) - S^T u_h^k = - [l^k (x)  (u)^*)]_(x_l^k)^(x_r^k) = - cal(E) ((u)^*, (u)^*)^T\
J^k M (dif u_h^k)/(dif t) - S^T u_h^k = - cal(E) ((u)^*, (u)^*)^T\
(dif u_h^k)/(dif t) = 1/J^k ( M^(-1) S^T u_h^k - M^(-1) cal(E) (f^*_R, f^*_L)^T)\
(dif u_h^k)/(dif t) = 1/J^k ( D_w u_h^k - "LIFT" (f^*_R, f^*_L)^T)
$

Now let's try multiple initial conditions with different number of elemnets and polynomial order.
First we will test non-smooth initial conditions. On the left are the initial
conditions and on the right is the solution after one period.

// cone advection
#figure(
  grid(
    columns: 2,
    gutter: 1mm,
    image("img/cone_12_4_initial.png", width: 100%), image("img/cone_12_4_final.png", width: 100%),
  ),
  caption: [ Cone with $K = 12$, $N = 4$ ],
  numbering: none,
)

#figure(
  grid(
    columns: 2,
    gutter: 1mm,
    image("img/cone_4_12_initial.png", width: 100%), image("img/cone_4_12_final.png", width: 100%),
  ),
  caption: [ Cone with $K = 4$, $N = 12$ ],
  numbering: none,
)

#figure(
  grid(
    columns: 2,
    gutter: 1mm,
    image("img/cone_12_12_initial.png", width: 100%), image("img/cone_12_12_final.png", width: 100%),
  ),
  caption: [ Cone with $K = 12$, $N = 12$ ],
  numbering: none,
)

Now let us test a discontinuous initial condition in the form of a step function.

#figure(
  grid(
    columns: 2,
    gutter: 1mm,
    image("img/saw_12_4_initial.png", width: 100%), image("img/saw_12_4_final.png", width: 100%),
  ),
  caption: [ Step signal with $K = 12$, $N = 4$ ],
  numbering: none,
)

#figure(
  grid(
    columns: 2,
    gutter: 1mm,
    image("img/saw_4_12_initial.png", width: 100%), image("img/saw_4_12_final.png", width: 100%),
  ),
  caption: [ Step signal with $K = 4$, $N = 12$ ],
  numbering: none,
)

#figure(
  grid(
    columns: 2,
    gutter: 1mm,
    image("img/saw_12_12_initial.png", width: 100%), image("img/saw_12_12_final.png", width: 100%),
  ),
  caption: [ Step signal with $K = 12$, $N = 12$ ],
  numbering: none,
)

== Solution reconstruction
- gibs phenomenon
- pointwise error
- pade reconstruction
- other reconstructions

// == Filters
// When approximating a discontinuous function using a polynomial we can observe the well known Gibbs phenomenon.

// = Nonlinearity
//
// == Limiters
//
// == Example nonlinear problem
// burger's equation
//
// = Systems of equations
//
// == Euler's equation

#bibliography("sources.bib")
