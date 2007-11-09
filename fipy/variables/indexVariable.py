## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "indexVariable.py"
 #                                     created: 10/25/07 {5:16:20 PM}
 #                                 last update: 11/9/07 {2:45:33 PM}
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

def _IndexVariable(index):
    if isinstance(index, tuple) or isinstance(index, list):
        return _ListIndexVariable(index)
    else:
        return _IndexVariable_(index)
                                  
class _IndexVariable_(Variable):
    def __init__(self, index):
        Variable.__init__(self, value=index, cached=1)
        
        self.index = self._requireIndex(index)
        
    def _requireIndex(self, index):
        return self._requires(self._checkIfSlice(index))

    def _checkIfSlice(x):
        if isinstance(x, slice):
            return _SliceVariable(x)
        else:
            return x
    _checkIfSlice = staticmethod(_checkIfSlice)
            
    def _indexValue(x):
        if isinstance(x, Variable):
            return x.getValue()
        else:
            return x
    _indexValue = staticmethod(_indexValue)
        
    def _calcValue(self):
        return self._indexValue(self.index)
        
    def _setValue(self, value, unit=None, array=None):
        self.value = value

    def _isItemMasked(x):
        return (isinstance(x, numerix.MA.MaskedArray)
                or (isinstance(x, Variable) and x._isMasked()))
    _isItemMasked = staticmethod(_isItemMasked)
    
    def _isMasked(self):
        return self._isItemMasked(self.index)
            
    def _maskedSlice(idx, a):
        def _filled(aa):
            def __filled(x):
                if isinstance(x, numerix.MA.MaskedArray):
                    return x.filled(0)
                else:
                    return x
                
            if isinstance(aa, tuple):
                return tuple([__filled(x) for x in aa])
            else:
                return __filled(aa)

        def _indexmask(aa, shape):
            (indexShape, 
             broadcastshape, 
             arrayindex) = numerix._indexShapeElements(idx, a.shape)
             
            if isinstance(aa, tuple):
                mask = numerix.MA.nomask
                for x in aa:
                    mask = numerix.MA.mask_or(mask, numerix.MA.getmask(x))
            else:
                mask = numerix.MA.getmask(aa)

            if arrayindex is not None and mask is not numerix.MA.nomask:
                tmp = numerix.MA.make_mask_none(indexShape[:arrayindex] 
                                                + broadcastshape
                                                + indexShape[arrayindex:])
                                                
                tmp[:] = mask[((numerix.newaxis,) * len(indexShape[:arrayindex])
                               + numerix.index_exp[:]
                               + (numerix.newaxis,) * len(indexShape[arrayindex:]))]
                               
                mask = tmp
            else:
                mask = numerix.MA.nomask
                
            return mask
            
        sliced = a[_filled(idx)]
        mask = _indexmask(idx, a.shape)
        
        if mask is not numerix.MA.nomask:
            mask = numerix.MA.mask_or(numerix.MA.getmask(sliced), mask)
            sliced = numerix.MA.array(data=sliced, mask=mask)
        elif numerix.MA.getmask(sliced) is numerix.MA.nomask:
            sliced = _filled(sliced)
            
        return sliced
    _maskedSlice = staticmethod(_maskedSlice)


    def _getitemfrom(self, other):
        if self._isMasked():
            op = lambda a, b: self._maskedSlice(a, b)
        else:
            op = lambda a, b: b[a]
        
        if not self.getUnit().isDimensionless():
            raise IndexError, "Indices must be dimensionless"
            
        if not isinstance(other, Variable):
            from fipy.variables.constant import _Constant
            other = _Constant(value=other)

        opShape, operatorClass, mesh = self._shapeClassAndMesh(None, None, other)
        
        if opShape is None or operatorClass is None:
            return NotImplemented
            
        from fipy.variables import binaryOperatorVariable
        binOp = binaryOperatorVariable._BinaryOperatorVariable(operatorClass)
        
        from fipy.variables.meshVariable import _MeshVariable

        if issubclass(binOp, _MeshVariable):
            return binOp(op=op, 
                         var=[self, other], 
                         mesh=mesh,
                         opShape=opShape, 
                         canInline=other.getUnit().isDimensionless(), 
                         unit=other.getUnit())
        else:
            return binOp(op=op, 
                         var=[self, other], 
                         opShape=opShape, 
                         canInline=other.getUnit().isDimensionless(), 
                         unit=other.getUnit())


    def _meshOperatorClass(self, opShape):
        from fipy.variables.meshVariable import _MeshVariable
        
        if isinstance(self.index, _MeshVariable) and opShape[-1] == self.index.shape[-1]:
            return self.index._OperatorVariableClass(), self.index.getMesh()
        else:
            return None, None

    def _shapeClassAndMesh(self, opShape, operatorClass, other):
        """
        Determine the shape of the result, the base class of the result, and (if
        necessary) a modified form of `other` that is suitable for the
        operation.
        """
        mesh = None
        
        if opShape is None:
            opShape = numerix._indexShape(index=self.getValue(), arrayShape=other.shape)
            
        if operatorClass is None:
            operatorClass, mesh = self._meshOperatorClass(opShape)
            if operatorClass is not None:
                return (opShape, operatorClass, mesh)

        opShape, baseClass, ignore = other._shapeClassAndOther(opShape, operatorClass, None)
        if baseClass is None:
            return (opShape, Variable._OperatorVariableClass(other, Variable), mesh)
        else:
            return (opShape, other._OperatorVariableClass(baseClass), mesh)
        
    def _repr(self, index):
        if isinstance(index, _SliceVariable):
            if index.start is None:
                s = ""
            else:
                s = repr(index.start)
            s += ":"
            if index.stop is not None:
                s += repr(index.stop)
            if index.step is not None:
                s += ":" + repr(index.step)
            return s
        elif index is Ellipsis:
            return "..."
        elif index is None:
            return "newaxis"
        else:
            return repr(index)
        
    def __repr__(self):
        return self._repr(self.index)
        
class _ListIndexVariable(_IndexVariable_):
    def _requireIndex(self, index):
        return [self._requires(self._checkIfSlice(i)) for i in index]

    def _calcValue(self):
        return tuple([self._indexValue(i) for i in self.index])
        
    def _isMasked(self):
        return numerix.logical_or.reduce([self._isItemMasked(i) for i in self.index])
        
    def _meshOperatorClass(self, opShape):
        from fipy.variables.meshVariable import _MeshVariable
        
        for item in self.index:
            if isinstance(item, _MeshVariable) and opShape[-1] == item.shape[-1]:
                return item._OperatorVariableClass(), item.getMesh()
                
        return None, None
        
    def __repr__(self):
        return ",".join([self._repr(i) for i in self.index])



        
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

