#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "iterator.py"
 #                                    created: 11/10/03 {2:47:38 PM} 
 #                                last update: 10/19/04 {2:52:01 PM} 
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
__docformat__ = 'restructuredtext'

import sys

class ConvergenceError(ArithmeticError):
    def __init__(self, equations):
	self.equations = equations
	
    def __str__(self):
	s = '\n'
	for equation in self.equations:
	    s += str(equation) + ' has residual = ' + str(equation._getResidual()) + '\n'
	return s

class Iterator:
    """Generic equation iterator
    """
    
    def __init__(self,equations,timeStepDuration = None):
	"""
	:Parameters:
	    
	      - `equations`: list or tuple of equations to iterate over
	      - `timeStepDuration`: duration of each timestep (`Variable`)
	"""
        self.equations = equations
	self.timeStepDuration = timeStepDuration
	
    def sweep(self, dt):
	for equation in self.equations:
	    equation.solve(dt)
	converged = True
	for equation in self.equations:
	    converged = converged and equation._isConverged()
	return converged
	
## 	if converged:
## 	    break
## 	elif maxSweeps > 1 and self.sweepCallback is not None:
## 	    self.sweepCallback(self.equations)
## ## 		print '\n'
## ## 		for equation in self.equations:
## ## 		    print str(equation) + ' has residual = ' + str(equation._getResidual())
## ## 		print '\n'
## ## 		equation.getVar().viewer.plot()
## 
## ## 		if (sweep + 1) % 10 == 0:
## ## 		    sys.stdout.write('|')
## ## 		else:
## ## 		    sys.stdout.write('.')
## ## 		sys.stdout.flush()
## 	    sweeping = 1
	    

	
    def sweeps(self, maxSweeps = 1, dt = 1.):
	converged = 0
	sweeping = 0
	for sweep in range(maxSweeps):
	    converged = self.sweep(dt)
	    
	    if converged:
		break
		
	if maxSweeps > 1 and not converged:
	    raise ConvergenceError(self.equations)
	
    def advanceTimeStep(self):
	for equation in self.equations:
	    var = equation.getVar()
	    var.updateOld()
	    
    def timestep(self, maxSweeps = 1, dt = 1.):
	"""Iterate the solution.
	
	:Parameters:
	    
	  - `maxSweeps`: maximum number of sweeps to reach convergence
	"""
	
	self.advanceTimeStep()
	self.sweeps(maxSweeps, dt)

    def elapseTime(self, desiredTime, maxSweepsPerStep = 1):
	elapsedTime = 0.
	while elapsedTime < desiredTime:
	    self.timestep(maxSweeps = maxSweepsPerStep)
	    elapsedTime += self.timeStepDuration
