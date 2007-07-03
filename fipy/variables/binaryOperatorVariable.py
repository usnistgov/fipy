## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "binaryOperatorVariable.py"
 #                                     created: 5/16/07 {9:55:54 AM}
 #                                 last update: 5/16/07 {9:55:54 AM}
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

from fipy.tools import numerix

def _BinaryOperatorVariable(operatorClass=None):
    # declare a binary operator class with the desired base class
    class binOp(operatorClass):

        def _calcValuePy(self):
            from fipy.variables.variable import Variable
            if isinstance(self.var[1], Variable):
                val1 = self.var[1].getValue()
            else:
                if type(self.var[1]) is type(''):
                    self.var[1] = physicalField.PhysicalField(value=self.var[1])
                val1 = self.var[1]

            return self.op(self.var[0].getValue(), val1)

        def getUnit(self):
            if self.unit is None:
                try:
                    return self._extractUnit(self.op(self.var[0]._getUnitAsOne(), self.var[1]._getUnitAsOne()))
                except:
                    return self._extractUnit(self._calcValuePy())
            else:
                return self.unit
                
        def _getRepresentation(self, style="__repr__", argDict={}, id=id, freshen=False):
            self.id = id
            return "(" + operatorClass._getRepresentation(self, style=style, argDict=argDict, id=id, freshen=freshen) + ")"
            
    return binOp

