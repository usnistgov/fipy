#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "adaptiveIterator.py"
 #                                    created: 11/10/03 {2:47:38 PM} 
 #                                last update: 7/19/04 {5:08:16 PM} 
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
    def __init__(self,equations,timeStepDuration = None,viewers = ()):
	"""Arguments:
	    
	    'equations' -- list or tuple of equations to iterate over
	"""
	Iterator.__init__(self, equations, timeStepDuration)
	self.viewers = viewers

    def resetTimeStep(self):
	for equation in self.equations:
	    var = equation.getVar()
	    var.resetToOld()
	    
    def adjustTimeStep(self, dt):
	factor = min([equation.getFigureOfMerit() for equation in self.equations])
	if factor > 1.:
	    factor = factor**0.5
	dt *= factor
	
	return dt

    def timestep(self, maxSweeps = 1, dt = 1.):
	self.elapsedTime = 0.
	self.desiredTime = dt
	while self.elapsedTime < self.desiredTime:
	    print "t:", self.elapsedTime, "dt:", dt
	    try:
		self.advanceTimeStep()
		self.sweeps(maxSweeps, dt)
		self.elapsedTime += dt
	    except ConvergenceError, e:
		print "ConvergenceError"
		self.resetTimeStep()
	    except KeyboardInterrupt:
		print "KeyboardInterrupt"
		break
## 	    except Exception, e:
## 		raise e
## 	    except:
## 		print "what's the error?"
		
		
	    dt = self.adjustTimeStep(dt)
		
	    for viewer in self.viewers:
		viewer.plot()

