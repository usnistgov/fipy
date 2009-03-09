## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "unaryOperatorVariable.py"
 #
 # Author: Jonathan Guyer <guyer@nist.gov>
 # Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 # Author: James Warren   <jwarren@nist.gov>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
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
 # #############################################################################
 ##

def _UnaryOperatorVariable(operatorClass=None):
    class unOp(operatorClass):
        def _calcValuePy(self):
            return self.op(self.var[0].getValue())

        def getUnit(self):
            if self.unit is None:
                try:
                    return self._extractUnit(self.op(self.var[0]._getUnitAsOne()))
                except:
                    return self._extractUnit(self._calcValue())
            else:
                return self.unit
            
    return unOp
