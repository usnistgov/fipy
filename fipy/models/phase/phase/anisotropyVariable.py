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
from fipy.tools.inline import inline
from fipy.variables.cellVariable import CellVariable
from fipy.variables.faceVariable import FaceVariable
from fipy.variables.vectorFaceVariable import VectorFaceVariable
from fipy.models.phase.phase.addOverFacesVariable import AddOverFacesVariable

class FFVariable(FaceVariable):
    def __init__(self, parameters = None, halfAngle = None):
        FaceVariable.__init__(self, halfAngle.getMesh())
        self.halfAngle = self.requires(halfAngle)
        self.parameters = parameters

    def calcValue(self):
        inline.optionalInline(self._calcValueIn, self._calcValuePy)

    def _calcValueIn(self):
        inline.runInlineLoop1("""
            double zsq = halfAngle(i) * halfAngle(i);
            double b = (1. - zsq) / (1. + zsq);
            double db = -N * 2. * halfAngle(i) / (1 + zsq);
            value(i) = alphasq * c2 * (1. + c2 * b) * db;
        """,halfAngle = self.halfAngle.getNumericValue(),
            N = self.parameters['symmetry'],
            value = self.value.value,
            alphasq = self.parameters['alpha']**2,
            c2 = self.parameters['anisotropy'],
            ni = len(self.value.value))
            
    def _calcValuePy(self):
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
        dPhi = self.phase.getFaceGrad()[:,:]
##        self.value = dPhi[:,::-1] * Numeric.array((-1.,1.))
        self.value[:,0] = -dPhi[:,1]
        self.value[:,1] = dPhi[:,0]

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
