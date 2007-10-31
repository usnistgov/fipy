## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "indexVariable.py"
 #                                     created: 10/25/07 {5:16:20 PM}
 #                                 last update: 10/31/07 {12:54:06 PM}
 # Author: Jonathan Guyer
 # E-mail: <jguyer@his.com>
 #   mail: Alpha Cabal
 #    www: <http://alphatcl.sourceforge.net>
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
 # History
 # 
 # modified   by  rev reason
 # ---------- --- --- -----------
 # 2007-10-25 JEG 1.0 original
 # 
 # #############################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.variables.variable import Variable
from fipy.tools import numerix

class _IndexVariable(Variable):
    def __init__(self, index):
        Variable.__init__(self, value=index, cached=1)
        
        def _checkIfSlice(x):
            if isinstance(x, slice):
                return _SliceVariable(x)
            else:
                return x
                
        if isinstance(index, tuple) or isinstance(index, list):
            self.index = [self._requires(_checkIfSlice(i)) for i in index]
        else:
            self.index = self._requires(_checkIfSlice(index))

    def _calcValue(self):
        if isinstance(self.index, list):
            value = []
            for i in self.index:
                if isinstance(i, Variable):
                    value.append(i.getValue())
                else:
                    value.append(i)
                    
            value = tuple(value)
        else:
            if isinstance(self.index, Variable):
                value = self.index.getValue()
            else:
                value = self.index
                
        return value
        
    def _setValue(self, value, unit=None, array=None):
        self.value = value

    def _isMasked(self):
        def __isMasked(x):
            return (isinstance(x, numerix.MA.MaskedArray)
                    or (isinstance(x, Variable) and x._isMasked()))
    
        if isinstance(self.index, list):
            return numerix.logical_or.reduce([__isMasked(i) for i in self.index])
        else:
            return __isMasked(self.index)

class _SliceVariable(Variable):
    def __init__(self, s):
        Variable.__init__(self, value=s, cached=1)
        self.start = self._requires(s.start)
        self.stop = self._requires(s.stop)
        self.step = self._requires(s.step)
        
    def _calcValue(self):
        def _getValue(x):
            if isinstance(x, Variable):
                return x.getValue()
            else:
                return x
            
        return slice(_getValue(self.start), 
                     _getValue(self.stop), 
                     _getValue(self.step))

    def _setValue(self, value, unit=None, array=None):
        self.value = value

