#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "testSteadyStateDiffusion.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 12/24/03 {10:18:38 AM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
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
from viewers.grid2DGistViewer import Grid2DGistViewer
import Numeric

def getParameters():
    L = 1.5
    def func(cell):
        r = L / 4.
        c = (L / 2., L / 2.)
        x = cell.getCenter()
        return (x[0] - c[0])**2 + (x[1] - c[1])**2 < r**2
    
    return {
        'nx'           :  100,
        'ny'           :  100,
        'L'            :  L,
        'theta value'  :  2. * Numeric.pi / 3.,
        'theta func'   :     func,
        'theta func value' : -2. * Numeric.pi / 3.
        }
    
if __name__ == '__main__':
    localParameters = getParameters()
    globalParameters = input.getParameters(localParameters)
    
    it = globalParameters['it']
    steps = globalParameters['steps']
    var = globalParameters['var']

    it.timestep(steps)

    viewer = Grid2DGistViewer(var)

    viewer.plot()
    raw_input()
            
