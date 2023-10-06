r"""Solve a single-component, liquid-vapor, van der Waals system.

This example solves a single-component, liquid-vapor, van der Waals system as
described by Wheeler *et al.* :cite:`PhysRevE.82.051601`. The free energy for this
system takes the form,

.. math:: :label: eq:reactiveWetting:liquidVapor1D:freeEnergy

   f = - \frac{e \rho^2}{m^2} + \frac{R T}{m} \left( \ln \frac{\rho}{m - \bar{v} \rho} \right)

where :math:`\rho` is the density. This free energy supports a two phase
equilibrium with densities given by :math:`\rho^l` and :math:`\rho^v` in the
liquid and vapor phases, respectively. The densities are determined by solving
the following system of equations,

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:pressureEquilibrium

   P \left( \rho^l \right) = P \left( \rho^v \right)

and

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:chemicalPotentialEquilibrium

   \mu \left( \rho^l \right) = \mu \left( \rho^v \right)

where :math:`\mu` is the chemical potential,

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:chemicalPotential

   \mu = \frac{\partial f}{\partial \rho}

and :math:`P` is the pressure,

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:pressure

   P = \rho \mu - f

One choice of thermodynamic parameters that yields a relatively physical two phase system is

>>> molarWeight = 0.118
>>> ee = -0.455971
>>> gasConstant = 8.314
>>> temperature = 650.
>>> vbar = 1.3e-05

with equilibrium density values of

>>> liquidDensity = 7354.3402662299995
>>> vaporDensity = 82.855803327810008

The equilibrium densities are verified by substitution into Eqs.
:eq:`eq:reactiveWetting:liquidVapor1D:pressureEquilibrium` and
:eq:`eq:reactiveWetting:liquidVapor1D:chemicalPotentialEquilibrium`. Firstly,
Eqs. :eq:`eq:reactiveWetting:liquidVapor1D:freeEnergy`,
:eq:`eq:reactiveWetting:liquidVapor1D:chemicalPotential` and
:eq:`eq:reactiveWetting:liquidVapor1D:pressure` are defined as python functions,

>>> from fipy import CellVariable, Grid1D, TransientTerm, VanLeerConvectionTerm, DiffusionTerm, ImplicitSourceTerm, ConvectionTerm, CentralDifferenceConvectionTerm, Viewer
>>> from fipy.tools import numerix

>>> def f(rho):
...     return ee * rho**2 / molarWeight**2 + gasConstant * temperature * rho / molarWeight * \
...            numerix.log(rho / (molarWeight - vbar * rho))

>>> def mu(rho):
...     return 2 * ee * rho / molarWeight**2 + gasConstant * temperature / molarWeight * \
...            (numerix.log(rho / (molarWeight - vbar * rho)) + molarWeight / (molarWeight - vbar * rho))

>>> def P(rho):
...     return rho * mu(rho) - f(rho)

The equilibrium densities values are verified with

>>> print(numerix.allclose(mu(liquidDensity), mu(vaporDensity)))
True

and

>>> print(numerix.allclose(P(liquidDensity), P(vaporDensity)))
True

In order to derive governing equations, the free energy functional is defined.

.. math::

   F = \int \left[ f + \frac{\epsilon T}{2} \left( \partial_j \rho \right)^2 \right] dV

Using standard dissipation laws, we write the governing equations for mass and
momentum conservation,

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:mass

   \frac{\partial \rho}{\partial t} + \partial_j \left(\rho u_j \right) = 0

and

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:momentum

   \frac{\partial \left( \rho u_i\right) }{\partial t} + \partial_j\left( \rho
   u_iu_j \right) = \partial_j\left( \nu \left[ \partial_j u_i + \partial_i u_j
   \right] \right) - \rho \partial_i \mu^{NC}

where the non-classical potential, :math:`\mu^{NC}`, is given by,

.. math::
  :label: eq:reactiveWetting:liquidVapor1D:nonClassicalPotential

   \mu^{NC} = \frac{\delta F}{\delta \rho} = \mu - \epsilon T \partial_j^2 \rho

As usual, to proceed, we define a mesh

>>> Lx = 1e-6
>>> nx = 100
>>> dx = Lx / nx
>>> mesh = Grid1D(nx=nx, dx=dx)

and the independent variables.

>>> density = CellVariable(mesh=mesh, hasOld=True, name=r'$\rho$')
>>> velocity = CellVariable(mesh=mesh, hasOld=True, name=r'$u$')
>>> densityPrevious = density.copy()
>>> velocityPrevious = velocity.copy()

The system of equations is solved in a fully coupled manner using a block
matrix. Defining :math:`\mu^{NC}` as an independent variable makes it easier to
script the equations without using higher order terms.

>>> potentialNC = CellVariable(mesh=mesh, name=r'$\mu^{NC}$')

>>> epsilon = 1e-16
>>> freeEnergy = (f(density) + epsilon * temperature / 2 * density.grad.mag**2).cellVolumeAverage

In order to solve the equations numerically, an interpolation method is used to
prevent the velocity and density fields decoupling. The following velocity
correction equation (expressed in discretized form) prevents decoupling from
occurring,

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:correction

   u_{i,f}^c = \frac{A_f d_f}{\overline{d}_f} \left( \overline{ \rho \partial_i
   \mu^{NC} }_f - \overline{\rho}_f \partial_{i,f} \mu^{NC} \right)

where :math:`A_f` is the face area, :math:`d_f` is the distance between the
adjacent cell centers and :math:`\overline{a}_f` is the momentum conservation
equation's matrix diagonal. The overbar refers to an averaged value between the
two adjacent cells to the face. The notation :math:`\partial_{i,f}` refers to a
derivative evaluated directly at the face (not averaged). The variable
:math:`u_i^c` is used to modify the velocity used in
Eq. :eq:`eq:reactiveWetting:liquidVapor1D:mass` such that,

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:massCorrected

   \frac{\partial \rho}{\partial t} + \partial_j \left(\rho \left[u_j + u_i^c
   \right] \right) = 0

Equation :eq:`eq:reactiveWetting:liquidVapor1D:massCorrected` becomes

>>> matrixDiagonal = CellVariable(mesh=mesh, name=r'$a_f$', value=1e+20, hasOld=True)
>>> correctionCoeff = mesh._faceAreas * mesh._cellDistances / matrixDiagonal.faceValue
>>> massEqn = TransientTerm(var=density) \
...           + VanLeerConvectionTerm(coeff=velocity.faceValue + correctionCoeff \
...                                         * (density * potentialNC.grad).faceValue, \
...                                   var=density) \
...           - DiffusionTerm(coeff=correctionCoeff * density.faceValue**2, var=potentialNC)

where the first term on the LHS of
Eq. :eq:`eq:reactiveWetting:liquidVapor1D:correction` is calculated in an
explicit manner in the ``VanLeerConvectionTerm`` and the second term is
calculated implicitly as a ``DiffusionTerm`` with :math:`\mu^{NC}` as the
independent variable.

In order to write Eq. :eq:`eq:reactiveWetting:liquidVapor1D:momentum` as a
:term:`FiPy` expression, the last term is rewritten such that,

.. math::

   \rho \partial_i \mu^{NC} = \partial_i \left( \rho \mu^{NC} \right) - \mu^{NC} \partial_i \rho

which results in

>>> viscosity = 1e-3
>>> ConvectionTerm = CentralDifferenceConvectionTerm
>>> momentumEqn = TransientTerm(coeff=density, var=velocity) \
...               + ConvectionTerm(coeff=[[1]] * density.faceValue * velocity.faceValue, var=velocity) \
...               == DiffusionTerm(coeff=2 * viscosity, var=velocity) \
...               - ConvectionTerm(coeff=density.faceValue * [[1]], var=potentialNC) \
...               + ImplicitSourceTerm(coeff=density.grad[0], var=potentialNC)

The only required boundary condition eliminates flow in or out of the domain.

>>> velocity.constrain(0, mesh.exteriorFaces)

As previously stated, the :math:`\mu^{NC}` variable will be solved implicitly. To
do this the Eq. :eq:`eq:reactiveWetting:liquidVapor1D:nonClassicalPotential` is
linearized in :math:`\rho` such that

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:chemicalPotentialEquation

   \mu^{NC} = \mu^* + \left( \frac{\partial \mu}{\partial \rho} \right)^* \left( \rho - \rho^* \right) - \epsilon T \partial_j^2 \rho

The :math:`^*` superscript denotes the current held value. In :term:`FiPy`,
:math:`\frac{\partial \mu}{\partial \rho}` is written as,

>>> potentialDerivative = 2 * ee / molarWeight**2 + gasConstant * temperature * molarWeight / density / (molarWeight - vbar * density)**2

and :math:`\mu^*` is simply,

>>> potential = mu(density)

Eq. :eq:`eq:reactiveWetting:liquidVapor1D:chemicalPotentialEquation` can be
scripted as

>>> potentialNCEqn = ImplicitSourceTerm(coeff=1, var=potentialNC) \
...                  == potential \
...                  + ImplicitSourceTerm(coeff=potentialDerivative, var=density) \
...                  - potentialDerivative * density \
...                  - DiffusionTerm(coeff=epsilon * temperature, var=density)

Due to a quirk in :term:`FiPy`, the gradient of :math:`\mu^{NC}` needs to be
constrained on the boundary.  This is because ``ConvectionTerm``'s will
automatically assume a zero flux, which is not what we need in this case.

>>> potentialNC.faceGrad.constrain(value=[0], where=mesh.exteriorFaces)

All three equations are defined and an are combined together with

>>> coupledEqn = massEqn & momentumEqn & potentialNCEqn

The system will be solved as a phase separation problem with an initial density
close to the average density, but with some small amplitude noise. Under these
circumstances, the final condition should be two separate phases of roughly equal
volume. The initial condition for the density is defined by

>>> numerix.random.seed(2011)
>>> density[:] = (liquidDensity + vaporDensity) / 2 * \
...    (1  + 0.01 * (2 * numerix.random.random(mesh.numberOfCells) - 1))

Viewers are also defined.

>>> from fipy import input
>>> if __name__ == '__main__':
...     viewers = Viewer(density), Viewer(velocity), Viewer(potentialNC)
...     for viewer in viewers:
...         viewer.plot()
...     input('Arrange viewers, then press <return> to proceed...')
...     for viewer in viewers:
...         viewer.plot()

The following section defines the required control parameters. The `cfl`
parameter limits the size of the time step so that 
`dt = cfl * dx / max(velocity)`.

>>> cfl = 0.1
>>> tolerance = 1e-1
>>> dt = 1e-14
>>> timestep = 0
>>> relaxation = 0.5
>>> if __name__ == '__main__':
...     totalSteps = 1e10
... else:
...     totalSteps = 10

In the following time stepping scheme a time step is recalculated if the residual
increases between sweeps or the required tolerance is not attained within 20
sweeps. The major quirk in this scheme is the requirement of updating the
``matrixDiagonal`` using the entire coupled matrix. This could be achieved more
elegantly by calling ``cacheMatrix()`` only on the necessary part of the
equation. This currently doesn't work properly in :term:`FiPy`.

>>> while timestep < totalSteps:
... 
...     sweep = 0
...     dt *= 1.1
...     residual = 1.
...     initialResidual = None
... 
...     density.updateOld()
...     velocity.updateOld()
...     matrixDiagonal.updateOld()
... 
...     while residual > tolerance:
... 
...         densityPrevious[:] = density
...         velocityPrevious[:] = velocity
...         previousResidual = residual
... 
...         dt = min(dt, dx / max(abs(velocity)) * cfl)
... 
...         coupledEqn.cacheMatrix()
...         residual = coupledEqn.sweep(dt=dt)
... 
...         if initialResidual is None:
...             initialResidual = residual
... 
...         residual = residual / initialResidual
... 
...         if residual > previousResidual * 1.1 or sweep > 20:
...             density[:] = density.old
...             velocity[:] = velocity.old
...             matrixDiagonal[:] = matrixDiagonal.old
...             dt = dt / 10.
...             if __name__ == '__main__':
...                 print('Recalculate the time step')
...             timestep -= 1
...             break
...         else:
...             matrixDiagonal[:] = coupledEqn.matrix.takeDiagonal()[mesh.numberOfCells:2 * mesh.numberOfCells]
...             density[:] = relaxation * density + (1 - relaxation) * densityPrevious
...             velocity[:] = relaxation * velocity + (1 - relaxation) * velocityPrevious
... 
...         sweep += 1
... 
...     if __name__ == '__main__' and timestep % 10 == 0:
...         print('timestep: %e / %e, dt: %1.5e, free energy: %1.5e' % (timestep, totalSteps, dt, freeEnergy))
...         for viewer in viewers:
...             viewer.plot()
... 
...     timestep += 1

>>> from fipy import input
>>> if __name__ == '__main__':
...     input('finished')

>>> print(freeEnergy < 1.5e9)
True

.. .. bibmissing:: /refs.bib
    :sort:
"""
from __future__ import unicode_literals

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
