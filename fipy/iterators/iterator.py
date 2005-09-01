#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "iterator.py"
 #                                    created: 11/10/03 {2:47:38 PM} 
 #                                last update: 8/23/05 {11:34:25 AM} 
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
    """
    This simple iterator
    """
    def __init__(self, solve, preSolve = None, postSolve = None):
        """
        `preSolve`, `solve`, and `postSolve` are static functions to be
        called before each solution iteration, for the iteration, and after
        iterating, respectively.  
        
        `solve()` is the only function that is required and must be of the
        form::
            
            def solve(dt):
                eq1.solve(var = var1, dt = dt, ...)
                eq2.solve(var = var2, dt = dt, ...)
                eq3.solve(var = var3, dt = dt, ...)
                  :
                  :
                
                return residual

        where you should actually solve the equations, sweeping if desired. 
        
        `residual` is a positive real number that you can determine any way
        you like, although something of the form::
            
            residual = max(max(eq1.residual / eq1.maxError), max(eq2.residual / eq2.maxError), ...)
            
        is typical.  If `residual <= 1.`, the iterator will accelerate and
        proceed to the next iteration, and if `residual > 1.`, it will
        decelerate and reattempt the same iteration.
        
        `preSolve()` is where you will call `updateOld()` on your solution
        variables, if desired.
        
        `postSolve()` is where you would adjust the timestep, if desired.
        """
        self._solve = solve
        if preSolve is not None:
            self._preSolve = preSolve
        if postSolve is not None:
            self._postSolve = postSolve
            
        self.residual = 0.
            
    def _preSolve():
        pass
    _preSolve = staticmethod(_preSolve)
    
    def _doPreSolve(self):
        self._preSolve()
        
    def _postSolve(dtActual, dtTry):
        pass
    _postSolve = staticmethod(_postSolve)
    
    def _doPostSolve(self, dtActual):
        self._postSolve(dtActual = dtActual, dtTry = self.dtTry)
    
    def _doSolve(self, dt):
        return self._solve(dt = dt())

    def step(self, dtTotal = 1., dtTry  = None, dtActual = None):
        """
        Takes one timestep of duration `dtTotal`, taking iterations of duration
        `dtTry`.
        
        Returns the next `dtTry` to attempt.
        """
        
        def makeVariable(x):
            if x is None:
                x = dtTotal
            from fipy.variables.variable import Variable
            if not isinstance(x, Variable):
                x = Variable(value = x)
            
            return x

        self.dtTotal = dtTotal
        self.dtTry = makeVariable(dtTry)
        dtActual = makeVariable(dtActual)
        
        self.elapsed = 0.
        
##         print "%20s | %20s | %20s" % ("elapsed", "next dt", "residual")
        while 1:
            dtActual.setValue(min(self.dtTry(), self.dtTotal - self.elapsed))
            
            self._doPreSolve()
            self.residual = self._doSolve(dt = dtActual)

            self.elapsed += dtActual()

            self._doPostSolve(dtActual = dtActual)

            if self.elapsed >= dtTotal:
                return self.dtTry
                
##             print "%20g | %20g | %20g" % (self.elapsed, self.dtTry, self.residual)
            import sys
            sys.stdout.flush()
