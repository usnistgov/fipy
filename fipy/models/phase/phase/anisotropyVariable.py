#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "anisotropyVariable.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 1/20/04 {11:30:33 AM} { 2:35:45 PM}
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

from fivol.variables.cellVariable import CellVariable
from fivol.variables.faceVariable import FaceVariable
from fivol.variables.vectorFaceVariable import VectorFaceVariable
from fivol.examples.phase.phase.addOverFacesVariable import AddOverFacesVariable

class FFVariable(FaceVariable):
    def __init__(self, parameters = None, halfAngle = None):
        FaceVariable.__init__(self, halfAngle.getMesh())
        self.halfAngle = self.requires(halfAngle)
        self.parameters = parameters
        
    def calcValue(self):
        alpha = self.parameters['alpha']
	N = self.parameters['symmetry']
	c2 = self.parameters['anisotropy']
        
        zsq = self.halfAngle[:] * self.halfAngle[:]
	b = (1. - zsq) / (1. + zsq)
	db = -N * 2 * self.halfAngle[:] / (1 + zsq)
        self.value = alpha**2 * c2 * (1. + c2 * b) * db

class DPhiReverse(VectorFaceVariable):
    def __init__(self, phase):
        VectorFaceVariable.__init__(self, phase.getMesh())
        self.phase = self.requires(phase)
        
    def calcValue(self):
        dphi = self.phase.getFaceGrad()[:,:]
        self.value = dphi[:,::-1] * Numeric.array((-1.,1))
    
class AnisotropyVariable(CellVariable):
    def __init__(self, parameters = None, phase = None, halfAngle = None):
        CellVariable.__init__(self, phase.getMesh())
        self.requires(phase)
        self.requires(halfAngle)
        ff = FFVariable(parameters = parameters, halfAngle = halfAngle)
        dPhiReverse = DPhiReverse(phase)
        
        self.AOF = AddOverFacesVariable(faceVariable = ff, faceGradient = dPhiReverse)
        
    def calcValue(self):
        self.value = self.AOF.getNumericValue()
##	alpha = self.parameters['alpha']
##	N = self.parameters['symmetry']
##	c2 = self.parameters['anisotropy']
##        dphi = self.phase.getFaceGrad()[:,:]
        
##        zsq = self.halfAngle[:] * self.halfAngle[:]
##	b = (1. - zsq) / (1. + zsq)
##	db = -N * 2 * self.halfAngle[:] / (1 + zsq)
##        ff = alpha**2 * c2 * (1. + c2 * b) * db

##        dphiReverse = dphi[:,::-1] * Numeric.array((-1.,1))

##        self.value = toolsTmp.addOverFaces(faceGradient = dphiReverse,
##                                        faceVariable = ff,
##                                        mesh = self.mesh,
##                                        NCells = len(self.value[:]))
        
##        contributions = Numeric.sum(self.mesh.getAreaProjections() * dphiReverse,1)

##        contributions = contributions * ff
        
##        NIntFac = len(self.mesh.getInteriorFaces())
##        NExtFac = len(self.mesh.getFaces()) - NIntFac
        
##        contributions = Numeric.concatenate((contributions[:NIntFac], Numeric.zeros(NExtFac,'d')))
##        ids = self.mesh.getCellFaceIDs()

##        contributions = Numeric.take(contributions, ids)

##        NCells = len(self.value[:])
##	NMaxFac = self.mesh.getMaxFacesPerCell()

##        contributions = Numeric.reshape(contributions,(NCells,-1))

##        orientations = Numeric.reshape(self.mesh.getCellFaceOrientations(),(NCells,-1))

##        self.value = Numeric.sum(orientations*contributions,1) / self.mesh.getCellVolumes()

##    def addOverFaces(grad, diffusion, mesh, Ncells):

##        contributions = Numeric.sum(mesh.getAreaProjections() * grad,1)

##        contributions = contributions * ff
        
##        NIntFac = len(mesh.getInteriorFaces())
##        NExtFac = len(mesh.getFaces()) - NIntFac
        
##        contributions = Numeric.concatenate((contributions[:NIntFac], Numeric.zeros(NExtFac,'d')))
##        ids = mesh.getCellFaceIDs()

##        contributions = Numeric.take(contributions, ids)

##	NMaxFac = mesh.getMaxFacesPerCell()

##        contributions = Numeric.reshape(contributions,(NCells,-1))

##        orientations = Numeric.reshape(mesh.getCellFaceOrientations(),(NCells,-1))

##        return Numeric.sum(orientations*contributions,1) / self.mesh.getCellVolumes()
        

        

