#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "adaptiveIterator.py"
 #                                    created: 11/10/03 {2:47:38 PM} 
 #                                last update: 4/2/04 {4:06:03 PM} 
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

"""Generic equation iterator
"""

from fipy.iterators.iterator import Iterator
from fipy.iterators.iterator import ConvergenceError


class AdaptiveIterator(Iterator):
    """Adaptive time-stepping equation iterator
    """
    
    def resetTimeStep(self):
	for equation in self.equations:
	    var = equation.getVar()
	    var.resetToOld()
	    
    def elapseTime(self, desiredTime, maxSweepsPerStep = 1):
	elapsedTime = 0.
	while elapsedTime < desiredTime:
	    print "t:", elapsedTime, "dt:", self.timeStepDuration
	    try:
		self.advanceTimeStep()
		self.sweeps(maxSweepsPerStep)
		elapsedTime += self.timeStepDuration.getValue()
	    except ConvergenceError, e:
		print "ConvergenceError"
		self.resetTimeStep()
	    except KeyboardInterrupt:
		print "KeyboardInterrupt"
		break
	    except Exception, e:
		print "Error:",e
	    except:
		print "what's the error?"
		
		
	    # adjust timestep
	    factors = [(equation.getSolutionTolerance()/equation.getResidual())**0.2 for equation in self.equations]
	    self.timeStepDuration.setValue(self.timeStepDuration.getValue() * min(factors))
		

