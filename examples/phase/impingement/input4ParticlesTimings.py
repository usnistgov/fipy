#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 1/16/04 {12:00:06 PM}
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
from fivol.examples.phase.examples.impingement.input4Particles import System4Particles
import Numeric
import time


def timings():
    
    resultsFile = open('results.txt', 'w')
    resultsFile.write('N     mesh build time      run time\n')
    resultsFile.close()


    for n in (10, 20, 40, 80, 160, 320, 640):
        for i in range(5):
            
            t0 = time.time()
            system = System4Particles(nx = n, ny = n, steps = 100, drivingForce = 10.)
            
            t1 = time.time()
            system.run()
            
            t2 = time.time()
            meshBuildTime = t1 - t0
            runTime =  t2 - t1
            
            resultsFile = open('results.txt', 'a')
            resultsFile.write('%i     %f       %f\n' % (n, meshBuildTime, runTime))
            resultsFile.close()

if __name__ == '__main__':
    timings()
