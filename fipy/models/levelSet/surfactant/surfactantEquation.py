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

        transientCoeff = Numeric.where(distanceVar < 0., 1e10, 0.)

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
        
    def calcValue(self):
##        id1, id2 = self.mesh.getAdjacentCellIDs()
##        epsilon = 1e-10

##        cellGrad = self.distanceVar.getGrad()
##        shape = Numeric.shape(cellGrad)
##        cellGradMag = Numeric.reshape(Numeric.repeat(cellGrad.getMag(), shape[1]), shape)

##        cellNormal = Numeric.where(cellGradMag > epsilon,
##                                   cellGrad / cellGradMag,
##                                   Numeric.zeros(shape))

##        faceNormal1 = array.take(cellNormal, id1)
##        faceNormal2 = array.take(cellNormal, id2)
##        var1 = array.take(self.distanceVar, id1)
##        var2 = array.take(self.distanceVar, id2)
##        shape = Numeric.shape(faceNormal1)
##        var1 = Numeric.reshape(Numeric.repeat(var1, shape[1]), shape)
##        var2 = Numeric.reshape(Numeric.repeat(var2, shape[1]), shape)

##        self.value = -Numeric.where(var2 > var1, faceNormal1, faceNormal2)
##        print "convectionCoeff",array.take(self.value, self.mesh.getCellFaceIDs())[3]
##        raw_input('convectionCoeff')
        
        ## interior faces
        faceGrad = self.distanceVar.getFaceGrad()
        faceGradMag = Numeric.array(faceGrad.getMag())
        faceGrad = Numeric.array(faceGrad)
        faceGradMag = Numeric.where(faceGradMag < 1e-10, 1e-10, faceGradMag)
        self.value = -faceGrad / faceGradMag[:, Numeric.NewAxis]

        ## exterior faces
        id1 = self.mesh.getAdjacentCellIDs()[0]
        cellGradMag = self.distanceVar.getGrad().getMag()
        cellGradMag = Numeric.where(cellGradMag < 1e-10, 1e-10, cellGradMag)
        cellNormal = self.distanceVar.getGrad() / cellGradMag[:, Numeric.NewAxis]
        faceCellNormals = Numeric.take(cellNormal, id1)
        mask = (self.mesh.getFaceCellIDs()[:,1]).mask()
        shape = Numeric.shape(faceCellNormals)
        mask = Numeric.reshape(Numeric.repeat(mask, shape[1]), shape)
        self.value = Numeric.where(mask, -faceCellNormals, self.value)

