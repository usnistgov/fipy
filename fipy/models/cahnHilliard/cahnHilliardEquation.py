#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "cahnHilliardEquation.py"
 #                                    created: 7/2/04 {10:39:23 AM} 
 #                                last update: 7/2/04 {11:51:58 AM} 
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

import Numeric

from fipy.equations.matrixEquation import MatrixEquation
from fipy.terms.transientTerm import TransientTerm
from fipy.terms.scSourceTerm import ScSourceTerm

class CahnHilliardEquation(MatrixEquation):

    def __init__(self,
                 var,
                 solver = 'default_solver',
                 boundaryConditions = (),
                 parameters = {}):
        
        mesh = var.getMesh()
	self.parameters = parameters
	self.var = var
        diffusionCoeff = self.parameters['diffusionCoeff']
        
        doubleWellDerivative = diffusionCoeff * asq * ( 1 - 2 * var**2 ) * (1 - var) * var
        laplacian = doubleWellDerivative.getLaplacian()
        
        scTerm = ScSourceTerm(
            sourceCoeff = laplacian,
	    mesh = var.getMesh())

        nthOrderTerm = NthOrderTerm(
            coeffs = (diffusionCoeff, self.parameters['epsilon']**2), 
            mesh = var.getMesh(),
            boundaryConditions = boundaryConditions)
            
	terms = (
	    TransientTerm(1., mesh),
	    nthOrderTerm,
            scTerm,
	)

	MatrixEquation.__init__(
            self,
            var,
            terms,
            solver)
