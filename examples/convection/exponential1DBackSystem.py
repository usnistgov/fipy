#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "exponential1DBackSystem.py"
 #                                    created: 12/16/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:06:46 PM} 
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
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

from fipy.terms.exponentialConvectionTerm import ExponentialConvectionTerm
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux

from convectionDiffusionSystem import ConvectionDiffusionSystem

class Exponential1DBackSystem(ConvectionDiffusionSystem):

    def __init__(self):
        self.L = 10.
        self.nx = 1000
        self.ny = 1
        self.diffCoeff = 1.
        self.convCoeff = (-10.,0.)
        self.sourceCoeff = 0.
        self.convectionScheme = ExponentialConvectionTerm
        ConvectionDiffusionSystem.__init__(self)

    def getBoundaryConditions(self):
        return (
            FixedValue(faces = self.mesh.getFacesLeft(),value = self.valueLeft),
            FixedValue(faces = self.mesh.getFacesRight(),value = self.valueRight),
            FixedFlux(faces = self.mesh.getFacesTop(),value = 0.),
            FixedFlux(faces = self.mesh.getFacesBottom(),value = 0.)
            )
        
if __name__ == '__main__':
    system = Exponential1DBackSystem()
    system.run()

            
            
