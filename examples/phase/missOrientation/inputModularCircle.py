#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "testSteadyStateDiffusion.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 1/16/04 {10:51:45 AM} 
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

from __future__ import nested_scopes
from input import PhaseSystem

import Numeric

class ModularCircleSystem(PhaseSystem):
    def __init__(self):        
        self.L = 1.5
        def func(cell):
            r = self.L / 4.
            c = (self.L / 2., self.L / 2.)
            x = cell.getCenter()
            return (x[0] - c[0])**2 + (x[1] - c[1])**2 < r**2
        self.func = func
        self.nx = 100
        self.ny = 100
        self.thetaValue = 2. * Numeric.pi / 3.
        self.thetaFuncValue = -2. * Numeric.pi / 3.
        PhaseSystem.__init__(self)
    
if __name__ == '__main__':
    system = ModularCircleSystem()
    system.run()

