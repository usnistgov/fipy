#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "sourceEquation.py"
 #                                    created: 10/7/04 {10:39:23 AM} 
 #                                last update: 10/7/04 {10:35:49 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
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
 
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.equations.matrixEquation import MatrixEquation
from fipy.terms.transientTerm import TransientTerm
from fipy.terms.scSourceTerm import ScSourceTerm
from fipy.terms.spSourceTerm import SpSourceTerm

class SourceEquation(MatrixEquation):
    """

    Equation consiting of a transient term, and source terms.
    
    """    
    def __init__(self,
                 var,
                 transientCoeff = 1.,
                 scCoeff = 0.,
                 spCoeff = 0.,
                 solver = LinearPCGSolver(tolerance = 1.e-15,
                                          steps = 1000),
                 boundaryConditions = (),
                 otherTerms = ()):
        
        mesh = var.getMesh()

	terms = (
	    TransientTerm(transientCoeff,mesh),
            ScSourceTerm(scCoeff, mesh),
            SpSourceTerm(spCoeff, mesh))

	MatrixEquation.__init__(
            self,
            var,
            terms,
            solver)

