#set text(font: "New Computer Modern", lang: "en")
#set heading(numbering: "1.")
#set document(
  title: [#align(center)[GD method for the 1D linear advection equation]],
  author: "Jachym Smid",
  date: auto
)


#title()

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
Now we have to define $u_h (x_(i plus 1/2))$ to do that we define suitable numerical flux.
$ u_h (x_(x plus 1/2)) = hat(u)_h (x_(x plus 1/2)) $
One such numerical flux could be the well known upwind ($a > 0$)
$ hat(u) (x_(x plus 1/2)) = u_h^- (x_(x plus 1/2)) = limits(lim)_(x arrow x_(j + 1/2)) u_h(x) $
The integral equation then becomes
$ integral_(I_j) (partial u_h)/(partial t) v dif x - integral_(I_j) (partial v)/(partial x) u_h dif x + [ a hat(u)_h (x_(j - 1/2)) v (x_(j - 1/2)) - a hat(u)_h (x_(j + 1/2)) v (x_(j + 1/2))] = 0 $
If we sum all the elements we get the following equation
$ sum_(j = 1)^N (integral_(I_j) (partial u_h)/(partial t) v dif x - integral_(I_j) (partial v)/(partial x) u_h dif x + [ a hat(u)_h (x_(j - 1/2)) v (x_(j - 1/2)) - a hat(u)_h (x_(j + 1/2)) v (x_(j + 1/2))]) = 0 $

== Basis functions

When choosing the basis functions for the space $V_h^k$ we have many options, for example
- Monomial basis
- Legendre polynomials
- Lagrange polynomials in Gauss points

Let us focus on the monomial basis.

=== Monomial basis

For easier computation we construct our basis on a so called refrence element $hat(I) = chevron.l -1, 1 chevron.r$
and then we map it onto each element $I_j$ using a affine transformation. Which is given by
$ x = x_j + h/2 xi, $
where $xi in hat(I)$ is a local variable on $hat(I)$.
The monomial basis then is
$ hat(phi)^(m) (xi) = xi^m, quad m = 1, 2, dots k $ 
For computation itself we need to evaluate integrals of a product of two basis functions on $I_j$, we do that by substitution onto the reference interval
$
integral_(I_j) phi_j^m phi_j^n dif x = integral_(-1)^1 hat(phi)^m hat(phi)^n h/2 dif xi =\
= h/2 integral_(-1)^1 xi^(m+n) dif xi = cases(h/(m+n+1) quad &: quad (m+n) percent 2 = 0,0 quad &: quad (m+n) percent 2 != 0)
$
These form the mass matrix $M.$

Next we will also need integrals involving a derivative of a basis function
$
integral_(I_j) phi_j^m 
$



== Time discretization

