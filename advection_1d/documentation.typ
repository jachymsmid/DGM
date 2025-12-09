#set text(font: "New Computer Modern", lang: "en")
#set heading(numbering: "1.")


#title("GDM method for the 1D linear advection equation")

= Problem overview
The linear advection equation is of the following form:
$ (partial u)/(partial t) + a (partial u)/(partial x) = 0 $
where $a$ is the advection speed and $u = u(x,t)$ is the unknown scalar function.
This equation is hyperbolic.

For the initial conditions problem we have to specify the initial conditions $u(x,t) bar.v_(t = 0) = u_0 (x)$, this can be solved analytically using the method of charcteristics.

= Discontinous Galerkin Method

== Spatial discretization
We split our domain (interval $I$) into elements $I_k = chevron.l x_(i), x_() chevron.r$

== Time discretization
