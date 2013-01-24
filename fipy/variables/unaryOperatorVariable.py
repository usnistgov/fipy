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

__docformat__ = 'restructuredtext'

__all__ = []

def _UnaryOperatorVariable(operatorClass=None):
    """
    Test BinOp pickling

    >>> from fipy import Grid1D, FaceVariable, CellVariable, dump, Variable
    >>> import os, sys
    >>> m = Grid1D()
    >>> vs = (CellVariable(mesh=m, value=2.), FaceVariable(mesh=m, value=3.), Variable(4))
    >>> tmp = []
    >>> for v in vs:
    ...     (f, n) = dump.write(-v)
    ...     tmp += [dump.read(n)]
    ...     if f is not None: 
    ...         os.close(f)
    ...         os.remove(n)
    >>> for v in tmp:
    ...     print v.__class__
    <class 'fipy.variables.cellVariable.CellVariable'>
    <class 'fipy.variables.faceVariable.FaceVariable'>
    <class 'fipy.variables.variable.Variable'>
    >>> print tmp[0].allclose(-2.)
    True
    >>> print tmp[1].allclose(-3.)
    True
    >>> print tmp[2].allclose(-4.)
    True
    """
    
    class unOp(operatorClass):
        def _calcValue_(self):
            return self.op(self.var[0].value)

        @property
        def unit(self):
            assert(hasattr(self, "_unit") == True)
            if self._unit is None:
                try:
                    return self._extractUnit(self.op(self.var[0]._unitAsOne))
                except:
                    return self._extractUnit(self._calcValue())
            else:
                return self._unit
            
    return unOp

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test()   
