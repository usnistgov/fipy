#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "iterator.py"
 #                                    created: 11/10/03 {2:47:38 PM} 
 #                                last update: 12/22/03 {5:00:37 PM} 
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
    
    def __init__(self,equations,maxSweeps = 1):
	"""Arguments:
	    
	    'equations' -- list or tuple of equations to iterate over
	"""
        self.equations = equations
	self.maxSweeps = maxSweeps
	
    def iterate(self,steps):
	"""Iterate the solution.
	
	Arguments:
	    
	    'steps' -- number of iteration time steps
	    
	    'timeStep' -- duration of each time step
	
	The 'preSolve()' method of each equation is called to do any
	necessary initialization, the equations are 'solve()' ed, then the
	'postSolve()' methods are called to do any necessary cleanup.
	"""
	for step in range(steps):
            for equation in self.equations:
                var = equation.getVar()
                var.updateOld()
	    for sweep in range(self.maxSweeps):
		for equation in self.equations:
		    equation.solve()
		converged = 1 # Because Andy is too lazy to update to a Python written since the Eisenhower administration
		for equation in self.equations:
		    converged = converged and equation.isConverged()
		if converged:
		    break

