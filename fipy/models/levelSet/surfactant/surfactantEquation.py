#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "SurfactantEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 9/3/04 {10:40:02 PM} 
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

from fipy.equations.matrixEquation import MatrixEquation
from fipy.terms.transientTerm import TransientTerm
from fipy.terms.upwindConvectionTerm import UpwindConvectionTerm
from fipy.terms.explicitUpwindConvectionTerm import ExplicitUpwindConvectionTerm
from convectionCoeff import ConvectionCoeff
from fipy.solvers.linearLUSolver import LinearLUSolver
from fipy.boundaryConditions.fixedValue import FixedValue

class SurfactantEquation(MatrixEquation):
    """

    A `SurfactantEquation` aims to evolve a surfactant on an interface
    defined by the zero level set of the `distanceVar`. The method
    should completely conserve the total coverage of surfactant.  The
    surfactant is only in the cells immediately in front of the
    advancing interface. The method only works for a positive velocity
    as it stands.
    
    """
    
    def __init__(self,
                 var,
                 distanceVar,
                 solver = LinearLUSolver(tolerance = 1e-15),
                 boundaryConditions = None):

	mesh = var.getMesh()

        if boundaryConditions == None:
            boundaryConditions = (FixedValue(mesh.getExteriorFaces(), 0),)

	transientTerm = TransientTerm(1. ,mesh)

        convectionTerm = ExplicitUpwindConvectionTerm(ConvectionCoeff(distanceVar), mesh, boundaryConditions)

	terms = (
	    transientTerm,
	    convectionTerm
            )
	    
	MatrixEquation.__init__(
            self,
            var,
            terms,
            solver)

    def solve(self, dt = 1.):
        MatrixEquation.solve(self, dt = 1.)

