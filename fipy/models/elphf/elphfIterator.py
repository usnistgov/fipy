#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "iterator.py"
 #                                    created: 11/10/03 {2:47:38 PM} 
 #                                last update: 9/3/04 {10:35:38 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
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

import Numeric

from fipy.iterators.iterator import Iterator
from fipy.iterators.adaptiveIterator import AdaptiveIterator

## import sys

class ElPhFIterator(Iterator):
    def __init__(self,equations,timeStepDuration = None,viewers = ()):
	"""Arguments:
	    
	    'equations' -- list or tuple of equations to iterate over
	"""
	Iterator.__init__(self, equations, timeStepDuration)
	self.viewers = viewers
	
    def sweep(self):
	converged = Iterator.sweep(self)
	
## 	if not converged:
	print '\n'
	for equation in self.equations:
	    residual = equation.getResidual()
	    print str(equation) + ' has residual = ' \
	    + str(Numeric.sqrt(Numeric.dot(residual,residual))/len(residual)) \
	    + ' with relaxation = ' \
	    + str(Numeric.sqrt(Numeric.dot(equation.relaxation,equation.relaxation))/len(equation.relaxation))
	print '\n'
	for viewer in self.viewers:
	    viewer.plot()

## 		if (sweep + 1) % 10 == 0:
## 		    sys.stdout.write('|')
## 		else:
## 		    sys.stdout.write('.')
## 		sys.stdout.flush()
	    
	return converged
	
