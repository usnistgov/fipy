#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "input4ParticlesProfile.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 4/2/04 {4:00:45 PM}
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
from fipy.examples.phase.examples.impingement.input4Particles import System4Particles
import Numeric
from fipy.profiler.profiler import calibrate_profiler
from fipy.profiler.profiler import Profiler

import sys

class System4ParticlesProfile(System4Particles):

    def run(self):
        fudge = calibrate_profiler(10000)
        self.it.timestep(1)
        profile = Profiler('profile.txt', fudge=fudge)
        for i in range(5):
            self.it.timestep(1)
            print "timestep: ",i        

        profile.stop()

	print 'finished run'
    
if __name__ == '__main__':
    system = System4ParticlesProfile(nx = 50, ny = 50)
    system.run()
    
    raw_input()

