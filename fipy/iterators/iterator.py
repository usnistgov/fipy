#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "iterator.py"
 #                                    created: 11/10/03 {2:47:38 PM} 
 #                                last update: 12/29/03 {2:39:31 PM} 
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

class Iterator:
    """Generic equation iterator
    """
    
    def __init__(self,equations):
	"""Arguments:
	    
	    'equations' -- list or tuple of equations to iterate over
	"""
        self.equations = equations
	
    def sweep(self, maxSweeps = 1):
	converged = 0
	for sweep in range(maxSweeps):
	    for equation in self.equations:
		equation.solve()
	    converged = 1 # Because Andy is too lazy to update to a Python written since the Eisenhower administration
	    for equation in self.equations:
		converged = converged and equation.isConverged()
	    if converged:
		break
		
	return converged
	
    def advanceTimeStep(self):
	for equation in self.equations:
	    var = equation.getVar()
	    var.updateOld()

    def timestep(self, steps = 1, maxSweeps = 1):
	"""Iterate the solution.
	
	Arguments:
	    
	    'steps' -- number of iteration time steps
	    
	    'timeStep' -- duration of each time step
	"""
	for step in range(steps):
	    self.advanceTimeStep()
	    self.sweep(maxSweeps)
