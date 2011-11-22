#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "adsorbingSurfactantEquation.py"
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
 # protection and is in the public domain.  PFM is an experimental
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

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.variables.cellVariable import CellVariable
from fipy.models.levelSet.surfactant.surfactantEquation import SurfactantEquation
from fipy.terms.implicitSourceTerm import ImplicitSourceTerm
from fipy.solvers import DefaultAsymmetricSolver, LinearPCGSolver

__all__ = ["AdsorbingSurfactantEquation"]

class _AdsorptionCoeff(CellVariable):
    def __init__(self, distanceVar, bulkVar, rateConstant):
        CellVariable.__init__(self, mesh = distanceVar.mesh)

        self.distanceVar = self._requires(distanceVar)
        self.bulkVar = self._requires(bulkVar)
        self.rateConstant = rateConstant
        self.dt = 0

        
    def _calcValue(self):
        return numerix.array(self.dt * self.bulkVar
                             * self.rateConstant * self._multiplier())

    def _updateDt(self, dt):
        self.dt = dt
        self._markStale()

class _AdsorptionCoeffInterfaceFlag(_AdsorptionCoeff):
    def _multiplier(self):
        return self.distanceVar._cellInterfaceFlag
    
class _AdsorptionCoeffAreaOverVolume(_AdsorptionCoeff):
    def _multiplier(self):
        return self.distanceVar.cellInterfaceAreas / self.mesh.cellVolumes

class _MaxCoeff(CellVariable):
    def __init__(self, distanceVar, vars = ()):
        CellVariable.__init__(self, mesh = distanceVar.mesh)
        self.vars = vars
        for var in self.vars:
            self._requires(var)
        self.distanceVar = self._requires(distanceVar)

    def _calcMax(self):
        total = 0
        for var in self.vars:
            total += numerix.array(var.interfaceVar)
        return numerix.array(total > 1) * self.distanceVar._cellInterfaceFlag

class _SpMaxCoeff(_MaxCoeff):
    def _calcValue(self):
        return 1e20 * self._calcMax()

class _ScMaxCoeff(_MaxCoeff):
    def _calcValue(self):
        val = self.distanceVar.cellInterfaceAreas / self.mesh.cellVolumes
        for var in self.vars[1:]:
            val -= self.distanceVar._cellInterfaceFlag * numerix.array(var)

        return 1e20 * self._calcMax() * numerix.where(val < 0, 0, val)

class AdsorbingSurfactantEquation(SurfactantEquation):
    r"""

    The `AdsorbingSurfactantEquation` object solves the
    `SurfactantEquation` but with an adsorbing species from some bulk
    value. The equation that describes the surfactant adsorbing is
    given by,

    .. math::

       \dot{\theta} = J v \theta + k c (1 - \theta - \theta_{\text{other}}) - \theta c_{\text{other}} k_{\text{other}} - k^- \theta

    where :math:`\theta`, :math:`J`, :math:`v`, :math:`k`, :math:`c`,
    :math:`k^-` and :math:`n` represent the surfactant coverage, the curvature,
    the interface normal velocity, the adsorption rate, the concentration in the
    bulk at the interface, the consumption rate and an exponent of consumption,
    respectively. The :math:`\text{other}` subscript refers to another
    surfactant with greater surface affinity.

    The terms on the RHS of the above equation represent conservation of
    surfactant on a non-uniform surface, Langmuir adsorption, removal of
    surfactant due to adsorption of the other surfactant onto non-vacant sites
    and consumption of the surfactant respectively. The adsorption term is added
    to the source by setting :math:` S_c = k c (1 - \theta_{\text{other}})` and
    :math:`S_p = -k c`. The other terms are added to the source in a similar
    way.

    The following is a test case:

    >>> from fipy.models.levelSet.distanceFunction.distanceVariable \
    ...     import DistanceVariable
    >>> from fipy.models.levelSet.surfactant.surfactantVariable \
    ...     import SurfactantVariable
    >>> from fipy.meshes import Grid2D
    >>> dx = .5
    >>> dy = 2.3
    >>> dt = 0.25
    >>> k = 0.56
    >>> initialValue = 0.1
    >>> c = 0.2
    
    >>> from fipy.meshes import Grid2D
    >>> mesh = Grid2D(dx = dx, dy = dy, nx = 5, ny = 1)
    >>> distanceVar = DistanceVariable(mesh = mesh, 
    ...                                value = (-dx*3/2, -dx/2, dx/2, 
    ...                                          3*dx/2,  5*dx/2),
    ...                                hasOld = 1)
    >>> surfactantVar = SurfactantVariable(value = (0, 0, initialValue, 0 ,0), 
    ...                                    distanceVar = distanceVar)
    >>> bulkVar = CellVariable(mesh = mesh, value = (c , c, c, c, c))
    >>> eqn = AdsorbingSurfactantEquation(surfactantVar = surfactantVar,
    ...                                   distanceVar = distanceVar,
    ...                                   bulkVar = bulkVar,
    ...                                   rateConstant = k)
    >>> eqn.solve(surfactantVar, dt = dt)
    >>> answer = (initialValue + dt * k * c) / (1 + dt * k * c)
    >>> print numerix.allclose(surfactantVar.interfaceVar, 
    ...                  numerix.array((0, 0, answer, 0, 0)))
    1

    The following test case is for two surfactant variables. One has more
    surface affinity than the other.

    >>> from fipy.models.levelSet.distanceFunction.distanceVariable \
    ...     import DistanceVariable
    >>> from fipy.models.levelSet.surfactant.surfactantVariable \
    ...     import SurfactantVariable
    >>> from fipy.meshes import Grid2D
    >>> dx = 0.5
    >>> dy = 2.73
    >>> dt = 0.001
    >>> k0 = 1.
    >>> k1 = 10.
    >>> theta0 = 0.
    >>> theta1 = 0.
    >>> c0 = 1.
    >>> c1 = 1.
    >>> totalSteps = 10
    >>> mesh = Grid2D(dx = dx, dy = dy, nx = 5, ny = 1)
    >>> distanceVar = DistanceVariable(mesh = mesh, 
    ...                                value = dx * (numerix.arange(5) - 1.5),
    ...                                hasOld = 1)
    >>> var0 = SurfactantVariable(value = (0, 0, theta0, 0 ,0), 
    ...                           distanceVar = distanceVar)
    >>> var1 = SurfactantVariable(value = (0, 0, theta1, 0 ,0), 
    ...                           distanceVar = distanceVar)
    >>> bulkVar0 = CellVariable(mesh = mesh, value = (c0, c0, c0, c0, c0))
    >>> bulkVar1 = CellVariable(mesh = mesh, value = (c1, c1, c1, c1, c1))

    >>> eqn0 = AdsorbingSurfactantEquation(surfactantVar = var0,
    ...                                    distanceVar = distanceVar,
    ...                                    bulkVar = bulkVar0,
    ...                                    rateConstant = k0)

    >>> eqn1 = AdsorbingSurfactantEquation(surfactantVar = var1,
    ...                                    distanceVar = distanceVar,
    ...                                    bulkVar = bulkVar1,
    ...                                    rateConstant = k1,
    ...                                    otherVar = var0,
    ...                                    otherBulkVar = bulkVar0,
    ...                                    otherRateConstant = k0)

    >>> for step in range(totalSteps):
    ...     eqn0.solve(var0, dt = dt)
    ...     eqn1.solve(var1, dt = dt)
    >>> answer0 = 1 - numerix.exp(-k0 * c0 * dt * totalSteps)
    >>> answer1 = (1 - numerix.exp(-k1 * c1 * dt * totalSteps)) * (1 - answer0)
    >>> print numerix.allclose(var0.interfaceVar, 
    ...                  numerix.array((0, 0, answer0, 0, 0)), rtol = 1e-2)
    1
    >>> print numerix.allclose(var1.interfaceVar, 
    ...                  numerix.array((0, 0, answer1, 0, 0)), rtol = 1e-2)
    1
    >>> dt = 0.1
    >>> for step in range(10):
    ...     eqn0.solve(var0, dt = dt)
    ...     eqn1.solve(var1, dt = dt)

    >>> x, y = mesh.cellCenters
    >>> check = var0.interfaceVar + var1.interfaceVar
    >>> answer = CellVariable(mesh=mesh, value=check)
    >>> answer[x==1.25] = 1.
    >>> print check.allequal(answer)
    True

    The following test case is to fix a bug where setting the adosrbtion
    coefficient to zero leads to the solver not converging and an eventual
    failure.

    >>> var0 = SurfactantVariable(value = (0, 0, theta0, 0 ,0), 
    ...                           distanceVar = distanceVar)
    >>> bulkVar0 = CellVariable(mesh = mesh, value = (c0, c0, c0, c0, c0))

    >>> eqn0 = AdsorbingSurfactantEquation(surfactantVar = var0,
    ...                                    distanceVar = distanceVar,
    ...                                    bulkVar = bulkVar0,
    ...                                    rateConstant = 0)

    >>> eqn0.solve(var0, dt = dt)
    >>> eqn0.solve(var0, dt = dt)
    >>> answer = CellVariable(mesh=mesh, value=var0.interfaceVar)
    >>> answer[x==1.25] = 0.
    
    >>> print var0.interfaceVar.allclose(answer)
    True

    The following test case is to fix a bug that allows the accelerator to
    become negative.

    >>> nx = 5
    >>> ny = 5
    >>> dx = 1.
    >>> dy = 1.
    >>> mesh = Grid2D(dx=dx, dy=dy, nx = nx, ny = ny)
    >>> x, y = mesh.cellCenters

    >>> disVar = DistanceVariable(mesh=mesh, value=1., hasOld=True)
    >>> disVar[y < dy] = -1
    >>> disVar[x < dx] = -1
    >>> disVar.calcDistanceFunction()

    >>> levVar = SurfactantVariable(value = 0.5, distanceVar = disVar)
    >>> accVar = SurfactantVariable(value = 0.5, distanceVar = disVar)

    >>> levEq = AdsorbingSurfactantEquation(levVar,
    ...                                     distanceVar = disVar,
    ...                                     bulkVar = 0,
    ...                                     rateConstant = 0)

    >>> accEq = AdsorbingSurfactantEquation(accVar,
    ...                                     distanceVar = disVar,
    ...                                     bulkVar = 0,
    ...                                     rateConstant = 0,
    ...                                     otherVar = levVar,
    ...                                     otherBulkVar = 0,
    ...                                     otherRateConstant = 0)

    >>> extVar = CellVariable(mesh = mesh, value = accVar.interfaceVar)

    >>> from fipy.models.levelSet.advection.higherOrderAdvectionEquation \
    ...     import buildHigherOrderAdvectionEquation
    >>> advEq = buildHigherOrderAdvectionEquation(advectionCoeff = extVar)

    >>> dt = 0.1

    >>> for i in range(50):
    ...     disVar.calcDistanceFunction()
    ...     extVar.value = (numerix.array(accVar.interfaceVar))
    ...     disVar.extendVariable(extVar)
    ...     disVar.updateOld()
    ...     advEq.solve(disVar, dt = dt)
    ...     levEq.solve(levVar, dt = dt)
    ...     accEq.solve(accVar, dt = dt)

    >>> print (accVar >= -1e-10).all()
    True
    """
    def __init__(self,
                 surfactantVar = None,
                 distanceVar = None,
                 bulkVar = None,
                 rateConstant = None,
                 otherVar = None,
                 otherBulkVar = None,
                 otherRateConstant = None,
                 consumptionCoeff = None):
        """
        Create a `AdsorbingSurfactantEquation` object.

        :Parameters:
          - `surfactantVar`: The `SurfactantVariable` to be solved for.
          - `distanceVar`: The `DistanceVariable` that marks the interface.
          - `bulkVar`: The value of the `surfactantVar` in the bulk.
          - `rateConstant`: The adsorption rate of the `surfactantVar`.
          - `otherVar`: Another `SurfactantVariable` with more surface affinity.
          - `otherBulkVar`: The value of the `otherVar` in the bulk.
          - `otherRateConstant`: The adsorption rate of the `otherVar`.
          - `consumptionCoeff`: The rate that the `surfactantVar` is consumed during deposition.
                             
        """

          

        SurfactantEquation.__init__(self, distanceVar = distanceVar)

        spCoeff = _AdsorptionCoeffInterfaceFlag(distanceVar, bulkVar, rateConstant)
        scCoeff = _AdsorptionCoeffAreaOverVolume(distanceVar, bulkVar, rateConstant)

        self.eq += ImplicitSourceTerm(spCoeff) - scCoeff

        self.coeffs = (scCoeff, spCoeff)

        if otherVar is not None:
            otherSpCoeff = _AdsorptionCoeffInterfaceFlag(distanceVar, otherBulkVar, otherRateConstant)
            otherScCoeff = _AdsorptionCoeffAreaOverVolume(distanceVar, -bulkVar * otherVar.interfaceVar, rateConstant)

            self.eq += ImplicitSourceTerm(otherSpCoeff) - otherScCoeff

            self.coeffs += (otherScCoeff,)
            self.coeffs += (otherSpCoeff,)

            vars = (surfactantVar, otherVar)
        else:
            vars = (surfactantVar,)

        spMaxCoeff = _SpMaxCoeff(distanceVar, vars)
        scMaxCoeff = _ScMaxCoeff(distanceVar, vars)

        self.eq += ImplicitSourceTerm(spMaxCoeff) - scMaxCoeff - 1e-40

        if consumptionCoeff is not None:
            self.eq += ImplicitSourceTerm(consumptionCoeff)

    def solve(self, var, boundaryConditions=(), solver=None, dt=None):
        """
        Builds and solves the `AdsorbingSurfactantEquation`'s linear system once.
        	
        :Parameters:
           - `var`: A `SurfactantVariable` to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.
           
	"""

        for coeff in self.coeffs:
            coeff._updateDt(dt)
        if solver is None:
            import fipy.solvers.solver
            if fipy.solvers.solver == 'pyamg':
                from fipy.solvers.pyAMG.linearGeneralSolver import LinearGeneralSolver
                solver = LinearGeneralSolver(tolerance=1e-15, iterations=2000)
            else:
                solver = LinearPCGSolver()
            
        SurfactantEquation.solve(self, var, boundaryConditions=boundaryConditions, solver=solver, dt=dt)

    def sweep(self, var, solver=None, boundaryConditions=(), dt=None, underRelaxation=None, residualFn=None):
        r"""
        Builds and solves the `AdsorbingSurfactantEquation`'s linear
        system once. This method also recalculates and returns the
        residual as well as applying under-relaxation.

        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations. 
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.
           - `underRelaxation`: Usually a value between `0` and `1` or `None` in the case of no under-relaxation

	"""
        for coeff in self.coeffs:
            coeff._updateDt(dt)
        if solver is None:
            solver = DefaultAsymmetricSolver()
        return SurfactantEquation.sweep(self, var, solver=solver, boundaryConditions=boundaryConditions, dt=dt, underRelaxation=underRelaxation, residualFn=residualFn)

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test() 
