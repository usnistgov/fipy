#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "anisotropyVariable.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 1/13/04 {12:59:56 PM} { 2:35:45 PM}
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

from variables.cellVariable import CellVariable
import Numeric
import toolsTmp

class AnisotropyVariable(CellVariable):
    def __init__(self, parameters = None, phase = None, halfAngle = None):
        CellVariable.__init__(self, phase.getMesh())
	self.parameters = parameters
	self.phase = self.requires(phase)
        self.halfAngle = self.requires(halfAngle)

    def calcValue(self):
	alpha = self.parameters['alpha']
	N = self.parameters['symmetry']
	c2 = self.parameters['anisotropy']
        dphi = self.phase.getFaceGrad()[:,:]
        
        zsq = self.halfAngle[:] * self.halfAngle[:]
	b = (1. - zsq) / (1. + zsq)
	db = -N*2*self.halfAngle/(1+zsq)
        ff = alpha**2 * c2 * (1. + c2 * b) * db

        dphiReverse = dphi[:,::-1] * Numeric.array((-1.,1))

        self.value = tools.addOverFaces(faceGradient = dphiReverse,
                                        faceVariable = ff,
                                        mesh = self.mesh,
                                        NCells = len(self.value[:]))
        
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
        

        

