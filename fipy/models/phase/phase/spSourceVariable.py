#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "spSourceVariable.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 1/28/04 {4:20:21 PM} 
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

class SpSourceVariable(CellVariable):
    def __init__(self, theta = None, mPhi = None, phase = None, parameters = None):
        CellVariable.__init__(self, mesh = theta.getMesh())
        self.theta = self.requires(theta)
        self.phase = self.requires(phase)                                   
        self.parameters = parameters
        self.mPhi = self.requires(mPhi)
    
    def  calcValue(self):
        s = self.parameters['s']
	epsilon = self.parameters['epsilon']

        thetaMag = self.theta.getGrad().getMag()[:]
        
        spSourceCoeff = self.mPhi[:] * (self.phase[:] - (self.mPhi[:] < 0.))
        spSourceCoeff += (2*s + epsilon**2 * thetaMag) * thetaMag

	self.value = spSourceCoeff
