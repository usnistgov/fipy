#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "liquidVapor1D.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 # ###################################################################
 ##

r"""

This example solves a single-component, liquid-vapor, van der Waals system as
described by Wheeler et *al.* [PhysRevE.82.051601]_. The free energy for this
system takes the form,

.. math:: :label: eq:reactiveWetting:liquidVapor1D:freeEnergy
   
   f = - \frac{e \rho^2}{m^2} + \frac{R T}{m} \left( \ln \frac{\rho}{m - \bar{v} \rho} \right)

where :math:`\rho` is the density. This free energy supports a two phase in
equilibrium with densities given by :math:`\rho^l` and :math:`\rho^v` in the
liquid and vapor phases, respectively. These densities can be determined by
solving the following system of equations,

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

The equilibrium densities can be verified by substitution into Eqs.
:eq:`eq:reactiveWetting:liquidVapor1D:pressureEquilibrium` and
:eq:`eq:reactiveWetting:liquidVapor1D:chemicalPotentialEquilibrium`. Firstly,
Eqs. :eq:`eq:reactiveWetting:liquidVapor1D:freeEnergy`,
:eq:`eq:reactiveWetting:liquidVapor1D:chemicalPotential` and
:eq:`eq:reactiveWetting:liquidVapor1D:pressure` need to be defined as python
functions,

>>> from fipy import *

>>> def f(rho):
...     return ee * rho**2 / molarWeight**2 + gasConstant * temperature * rho / molarWeight * \
...            numerix.log(rho / (molarWeight - vbar * rho))

>>> def mu(rho):
...     return 2 * ee * rho / molarWeight**2 + gasConstant * temperature / molarWeight * \
...            (numerix.log(rho / (molarWeight - vbar * rho)) + molarWeight / (molarWeight - vbar * rho))

>>> def P(rho):
...     return rho * mu(rho) - f(rho)

The equilibrium densities values can then be verified.

>>> print numerix.allclose(mu(liquidDensity), mu(vaporDensity))
True

>>> print numerix.allclose(P(liquidDensity), P(vaporDensity))
True

In order to derive governing equations, the free energy functional must aslo be defined.

.. math::
   
   F = \int \left[ f + \frac{\epsilon T}{2} \left( \partial_j \rho \right)^2 \right] dV

Using standard dissipation laws, we can now write down the governing equations
for mass and momentum conservation,

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:continuity
   
   \frac{\partial \rho}{\partial t} + \partial_j \left(\rho u_j \right) = 0

and

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:momentum
   
   \frac{\partial \left( \rho u_i\right) }{\partial t} + \partial_j\left( \rho
   u_iu_j \right) = \partial_j\left( \nu \left[ \partial_j u_i + \partial_i u_j
   \right] \right) - \rho \partial_i \mu^{NC}

where the non-classical potential, :math:`\mu^{NC}` is given by,

.. math::
  :label: eq:reactiveWetting:liquidVapor1D:nonClassicalPotential

   \mu^{NC} = \frac{\delta F}{\delta \rho} = \mu - \epsilon T \partial_j^2 \rho
   
As usual, to proceed, we need to define a mesh

>>> Lx = 1e-6
>>> nx = 100
>>> dx = Lx / nx
>>> mesh = Grid1D(nx=nx, dx=dx)

and the independent variables.

>>> density = CellVariable(mesh=mesh, hasOld=True, name='$\rho$')
>>> velocity = CellVariable(mesh=mesh, hasOld=True, name='$u$')

The system of equations will be solved in a fully coupled manner using a block
matrix. Defining :math:`\mu^{NC}` as an independent variable makes it easier to
script the equations without using higher order terms.

>>> potentialNC = CellVariable(mesh=mesh, hasOld=True, name='$\mu^{NC}$')

In order to solve the equations numerically, an interpolation method must be used
to prevent the velocity and density fields decoupling. The following velocity
correction equation (expressed in discretized form) prevents decoupling from
occuring,

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
Eq. :eq:`eq:reactiveWetting:liquidVapor1D:continuity` such that,

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:continuityCorrected
   
   \frac{\partial \rho}{\partial t} + \partial_j \left(\rho \left[u_j + u_i^c
   \right] \right) = 0

Equation :eq:`eq:reactiveWetting:liquidVapor1D:continuityCorrected` can be
scripted in the form,

>>> matrixDiagonal = CellVariable(mesh=mesh, name='$a_f$')
>>> correctionCoeff = mesh._getFaceAreas() * mesh._getCellDistances() / matrixDiagonal.getFaceValue()
>>> densityEqn = TransientTerm(var=density) \
...              + VanLeerConvectionTerm(coeff=velocity.getFaceValue() + correctionCoeff \
...                                            * (density * potentialNC.getGrad()).getFaceValue(), \
...                                      var=density) \
...              - DiffusionTerm(coeff=correctionCoeff * faceDensity**2, var=potentialNC)

where the first term on the LHS of
Eq. :eq:`eq:reactiveWetting:liquidVapor1D:correction` is calculated in an
explicit manner in the ``VanLeerConvectionTerm`` and the second term is
calculated implicitly as a ``DiffusionTerm`` with :math:`\mu^{NC}` as the
independent variable.

In order to write Eq. :eq:`eq:reactiveWetting:liquidVapor1D:momentum` as a
:term:`FiPy` expression, the last term must be rewritten such that,

.. math::

   \rho \partial_i \mu^{NC} = \partial_i \left( \rho \mu^{NC} \right) - \mu^{NC} \partial_i \rho

which results in

>>> ConvectionTerm = CentralDifferenceConvectionTerm
>>> velocityEqn = TransientTerm(coeff=density, var=velocity) \
...               + ConvectionTerm(coeff=[[1]] * density.getFaceValue() * velocity.getFaceValue(), var=velocity) \
...               = DiffusionTerm(coeff=2 * viscosity, var=velocity) \
...               - ConvectionTerm(coeff=density.getFaceValue() * [[1]], var=potentialNC) \
...               + ImplicitSourceTerm(coeff=density.getGrad()[0], var=potentialNC)

The only required bnoundary condition eliminates flow in or out of the domain.

>>> velocity.constrain(0, mesh.getExteriorFaces())

As previously stated, the :math:`\mu^{NC}` variable will be solved implicitly. To
do this the Eq. eq:`eq:reactiveWetting:liquidVapor1D:nonClassicalPotential` needs
to be linearized in :math:`\rho` such that

.. math::
   :label: eq:reactiveWetting:liquidVapor1D:chemicalPotentialEquation

   \mu^{NC} = \mu^* + \left( \frac{\partial \mu}{\partial \rho} \right)^* \left( \rho - \rho^* \right) - \epsilon T \partial_j^2 \rho

The :math:`^*` superscript denotes the current held value. In :term:`FiPy`, can be written as,

>>> potentialDerivative = 2 * ee / molarWeight + gasConstant * temperature / density / (molarWeight - vbar * density)**2

and :math:`\mu^*` is simply,

>>> potential = mu(density)

Eq. :eq:`eq:reactiveWetting:liquidVapor1D:chemicalPotentialEquation` can then be
written as

>>> potententialNCEq = ImplicitSourceTerm(coeff=1, var=potentialNC) \
...                    == potential \
...                    + ImplicitSourceTerm(coeff=potentialDerivative, var=density) \
...                    - potentialDerivative * density \
...                    - DiffusionTerm(coeff=epsilon * temperature, var=density)

All three equations have now been defined and can now be combined together,

>>> coupledEq = densityEqn & velocityEq & potentialEq

The system will be solved as a phase separation problem with an initial density
close to the average density, but with some small amplitude noise. Under these
circumstances, the final condition should be two separate phases of roughly equal
volume. Define an initial condition for the density, such that

>>> density[:] = (liquidDensity + vaporDensity) / 2 * \
...    (1  + 0.01 * (2 * numerix.random.random(mesh.getNumberOfCells()) - 1))

>>> viewers = Viewer(density), Viewer(velocity), Viewer(potential)
>>> for viewer in viewers:
...     viewer.plot()

Some control parameters need to be defined. The ``cfl`` parameter limits the size
of the time step so that ``dt = cfl * dx / max(velocity)``. 

>>> cfl = 0.1
>>> tolerance = 1e-3
>>> globalTolerance = 1e-3
>>> dt = 1e-10
>>> globalResidual = 1.
>>> timestep = 0
>>> initialGlobalResidual = coupledEqn.justResidualVector(dt=1e+20)


>>> while globalResidual > globalTolerance: 
... 
...     residual = 1.
...     sweep = 0
...     dt *= 1.1
...     
...     density.updateOld()
...     velocity.updateOld()
...
...     globalResidual = coupledEqn.justResidualVector(dt=1e+20) / initialGlobalResidual
...     initialResiudal = coupledEqn.justResidualVector(dt=dt)
...
...     while residual > tolerance:
... 
...         dt = min(dt * 1.1, dx / max(abs(velocity)) * cfl)
... 
...         coupledEqn.cacheMatrix()
...         residual = coupledEqn.sweep(dt=dt) / initialResidual
...         ap[:] = coupledEqn.getMatrix().takeDiagonal()[mesh.getNumberOfCells(): 2 * mesh.getNumberOfCells()]
... 
...     for viewer in viewers:
...         viewer.plot()

"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
