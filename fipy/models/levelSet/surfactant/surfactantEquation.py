#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "SurfactantEquation.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 4/2/04 {4:00:26 PM} 
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

import MA
import Numeric

import fipy.tools.array as array

from fipy.equations.matrixEquation import MatrixEquation
from fipy.terms.transientTerm import TransientTerm
from fipy.terms.upwindConvectionTerm import UpwindConvectionTerm
from fipy.variables.vectorFaceVariable import VectorFaceVariable


class SurfactantEquation(MatrixEquation):
    def __init__(self,
                 var,
                 distanceVar,
                 solver='default_solver',
                 boundaryConditions=()):
		     
	mesh = var.getMesh()

        transientCoeff = (distanceVar < 0.) * 1e10

##        Numeric.where(distanceVar < 0., 1e10, 0.)

	transientTerm = TransientTerm(transientCoeff,mesh)
	convectionTerm = UpwindConvectionTerm(ConvectionCoeff(distanceVar), mesh, boundaryConditions)

	terms = (
	    transientTerm,
	    convectionTerm
            )
	    
	MatrixEquation.__init__(
            self,
            var,
            terms,
            solver)


class ConvectionCoeff(VectorFaceVariable):
    def __init__(self, distanceVar):
        VectorFaceVariable.__init__(self, distanceVar.getMesh(), name = 'surfactant convection')
        self.distanceVar = self.requires(distanceVar)

    def _calcValue(self):    
        faceGrad = self.distanceVar.getGrad().getArithmeticFaceValue()
        faceGradMag = Numeric.where(faceGrad.getMag() > 1e-10,
                                    faceGrad.getMag(),
                                    1e-10)
        faceGrad = Numeric.array(faceGrad)
        self.value = -faceGrad / faceGradMag[:,Numeric.NewAxis]

