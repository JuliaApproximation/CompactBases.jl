# Finite-differences overview

At the moment, there are three variations of finite-differences
implemented:

- [`FiniteDifferences`](@ref), standard, explicit, three-point
  stencil, uniform node distribution,
- [`StaggeredFiniteDifferences`](@ref), explicit, variationally
  derived three-point stencil, both uniform and non-uniform
  distributions. The nodes are staggered, such that in the uniform
  case, they follow the formula
  ```math
  r_j = (j-1/2)\rho
  ```
  where ``\rho`` is the step-size. This scheme is useful in polar and
  spherical coordinates. See

  - Schafer, K. J., Gaarde, M. B., Kulander, K. C., Sheehy, B., &
      DiMauro, L. F. (2000). Calculations of Strong Field Multiphoton
      Processes in Alkali Metal Atoms. AIP Conference Proceedings, 525(1),
      45–58. <http://dx.doi.org/10.1063/1.1291925>

  - Krause, J. L., & Schafer, K. J. (1999). Control of THz Emission from
      Stark Wave Packets. The Journal of Physical Chemistry A, 103(49),
      10118–10125. <http://dx.doi.org/10.1021/jp992144>

  - Koonin, S. E., & Meredith, D. C. (1990). Computational Physics,
      FORTRAN Version. Reading, Mass: Addison-Wesley.

  for the derivation.
- [`ImplicitFiniteDifferences`](@ref), also known as Compact
  Finite-Differences, where differential operators are approximated as
  ```math
  \vec{f}^{(o)}\approx\mat{M}_o^{-1}\Delta_o\vec{f}
  ```
  At present, only three-point stencils and uniform grids are supported.

  - Lele, S. K. (1992). Compact Finite Difference Schemes With
    Spectral-Like Resolution. Journal of Computational Physics, 103(1),
    16–42. <http://dx.doi.org/10.1016/0021-9991(92)90324-r>

  - Muller, H. G. (1999). An Efficient Propagation Scheme for the
    Time-Dependent Schrödinger equation in the Velocity Gauge. Laser
    Physics, 9(1), 138–148.

## Boundary conditions

At the moment, support for anything but Dirichlet0 boundary conditions
is lacking. Since the library was developed mostly with calculation
in atomic physics in mind, there is some special support for singular
Coulomb potentials built-in, which will remain there until the
functionality has been sufficiently generalized.
