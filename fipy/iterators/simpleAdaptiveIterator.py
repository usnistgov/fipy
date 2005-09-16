## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "rungeKutta5thOrderIterator.py"
 #                                    created: 8/12/05 {9:49:53 AM} 
 #                                last update: 9/16/05 {12:27:37 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  rungeKutta5thOrderIterator.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2005-08-12 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.iterators.iterator import Iterator, ConvergenceError

class SimpleAdaptiveIterator(Iterator):
    """
    Adaptive time-stepping equation iterator accelerates and decelerates as
    necessary to reach the total elapsed `dtTotal` in the fewest steps
    with acceptable error.
    """
    
    _SAFETY = 0.9
    _ERRCON = 1.89e-4
    
    def _doPostSolve(self, dtActual):
        if self.residual > self._ERRCON:
            self.dtTry.setValue(dtActual * self._SAFETY * self.residual**-0.25)
        else:
            self.dtTry.setValue(5. * dtActual)
            
        Iterator._doPostSolve(self, dtActual)
        
    def _doSolve(self, dt):
        while 1:
            residual = self._solve(dt = dt())
                
            if residual <= 1.:
                break     # step succeeded

            dt.setValue(dt * max(self._SAFETY * residual**-0.2, 0.1))
            
            if self.elapsed + dt() == self.elapsed:
                raise Exception, "step size underflow"

        print "dt:", dt

        return residual
