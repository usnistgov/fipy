#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 4/2/04 {4:06:31 PM}
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
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

from __future__ import nested_scopes
from examples.phase.impingement.input import ImpingementSystem

class System1D(ImpingementSystem):

    def initialConditions(self, Lx = None, Ly = None):

        cells = self.mesh.getCells()
        self.phase.setValue(1., cells)
        self.theta.setValue(1., cells)

        def getRightCells(cell, Lx = None):
            if cell.getCenter()[0] > Lx / 2.:
                return 1
            
        self.theta.setValue(0., self.mesh.getCells(getRightCells, Lx = Lx))
        
if __name__ == '__main__':
    system = System1D(nx = 40, ny = 1)
    system.run()
    raw_input()


