#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "SurfactantEquation.py"
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

__docformat__ = 'restructuredtext'

from fipy.terms.transientTerm import TransientTerm
from fipy.terms.explicitUpwindConvectionTerm import ExplicitUpwindConvectionTerm
from fipy.solvers import DefaultAsymmetricSolver

from fipy.models.levelSet.surfactant.convectionCoeff import _ConvectionCoeff

__all__ = ["SurfactantEquation"]

class SurfactantEquation:
    """

    A `SurfactantEquation` aims to evolve a surfactant on an interface
    defined by the zero level set of the `distanceVar`. The method
    should completely conserve the total coverage of surfactant.  The
    surfactant is only in the cells immediately in front of the
    advancing interface. The method only works for a positive velocity
    as it stands.
    
    """

    def __init__(self, distanceVar = None):
        """
        Creates a `SurfactantEquation` object.

        :Parameters:
          - `distanceVar`: The `DistanceVariable` that marks the interface.

        """
        transientTerm = TransientTerm(coeff = 1)

        convectionTerm = ExplicitUpwindConvectionTerm(_ConvectionCoeff(distanceVar))

        self.eq = transientTerm - convectionTerm

    def solve(self, var, boundaryConditions = (), solver=None, dt=None):
        """
        Builds and solves the `SurfactantEquation`'s linear system once.
                
        :Parameters:
           - `var`: A `SurfactantVariable` to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.

        """
        if type(boundaryConditions) not in (type(()), type([])):
            boundaryConditions = (boundaryConditions,)
        if solver is None:
            solver=DefaultAsymmetricSolver()

        var.constrain(0, var.mesh.exteriorFaces)

        self.eq.solve(var,
                      boundaryConditions=boundaryConditions,
                      solver = solver,
                      dt=1.)

    def sweep(self, var, solver=None, boundaryConditions=(), dt=None, underRelaxation=None, residualFn=None):
        r"""
        Builds and solves the `Term`'s linear system once. This method
        also recalculates and returns the residual as well as applying
        under-relaxation.

        :Parameters:

           - `var`: The variable to be solved for. Provides the initial condition, the old value and holds the solution on completion.
           - `solver`: The iterative solver to be used to solve the linear system of equations.
           - `boundaryConditions`: A tuple of boundaryConditions.
           - `dt`: The time step size.
           - `underRelaxation`: Usually a value between `0` and `1` or `None` in the case of no under-relaxation

	"""

        if type(boundaryConditions) not in (type(()), type([])):
            boundaryConditions = (boundaryConditions,)
        if solver is None:
            solver=DefaultAsymmetricSolver()
        
        var.constrain(0, var.mesh.exteriorFaces)

        return self.eq.sweep(var, solver=solver, boundaryConditions=boundaryConditions, underRelaxation=underRelaxation, residualFn=residualFn, dt=1.)
