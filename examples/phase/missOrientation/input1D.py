#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "input1D.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 1/16/04 {12:04:29 PM} 
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

import input
from input import PhaseSystem

class Phase1DSystem(PhaseSystem):
    def __init__(self):        
        self.L = 1.5
        self.func = lambda cell: cell.getCenter()[0] > self.L / 2.
        self.nx = 100
        self.ny = 1
        self.thetaValue = 1.
        self.thetaFuncValue = 0.
        PhaseSystem.__init__(self)
    
if __name__ == '__main__':
    system = Phase1DSystem()
    system.run()
