#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "input1DpoissonRightCharge.py"
 #                                    created: 1/15/04 {3:45:27 PM} 
 #                                last update: 5/5/04 {6:41:25 PM} 
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
 #  2004-01-15 JEG 1.0 original
 # ###################################################################
 ##

from fipy.tools.profiler.profiler import Profiler
from fipy.tools.profiler.profiler import calibrate_profiler

from fipy.viewers.grid2DGistViewer import Grid2DGistViewer

import input1Dpoisson

mesh = input1Dpoisson.mesh
it = input1Dpoisson.it
fields = input1Dpoisson.fields
parameters = input1Dpoisson.parameters
L = input1Dpoisson.L

setCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > L/2.)
fields['interstitials'][0].setValue(0.)
fields['interstitials'][0].setValue(1.,setCells)

if __name__ == '__main__':
    viewers = [Grid2DGistViewer(var = field) for field in fields['all']]

    for viewer in viewers:
	viewer.plot()
	
    raw_input()

    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)

    it.timestep(1)
    
    for viewer in viewers:
	viewer.plot()
	
    print fields['potential']
    print fields['potential'][0], fields['potential'][-1]
	
    # profile.stop()
	    
    raw_input()

