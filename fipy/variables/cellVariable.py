#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cellVariable.py"
 #                                    created: 12/9/03 {2:03:28 PM} 
 #                                last update: 1/3/07 {3:16:57 PM} 
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.variables.variable import Variable
from fipy.tools import numerix

        
class CellVariable(Variable):
    """
    Represents the field of values of a variable on a `Mesh`.
    
    A `CellVariable` can be ``pickled`` to persistent storage (disk) for later use:
        
        >>> from fipy.meshes.grid2D import Grid2D
        >>> mesh = Grid2D(dx = 1., dy = 1., nx = 10, ny = 10)
        
        >>> var = CellVariable(mesh = mesh, value = 1., hasOld = 1, name = 'test')
        >>> var.setValue(mesh.getCellCenters()[...,0] * mesh.getCellCenters()[...,1])

        >>> from fipy.tools import dump        
        >>> (f, filename) = dump.write(var, extension = '.gz')
        >>> unPickledVar = dump.read(filename, f)
        
        >>> print var.allclose(unPickledVar, atol = 1e-10, rtol = 1e-10)
        1
        
    """
    
    def __init__(self, mesh, name = '', value=0., unit = None, hasOld = 0):
        if value is None:
            array = None
        else:
            array =  numerix.zeros(self._getShapeFromMesh(mesh),'d')
#       array[:] = value
        
        Variable.__init__(self, mesh = mesh, name = name, value = value, unit = unit, array = array)

        if hasOld:
            self.old = self.copy()
        else:
            self.old = None
            
        self.arithmeticFaceValue = self.harmonicFaceValue = self.grad = self.faceGrad = self.volumeAverage = None
        
    def _getVariableClass(self):
        return CellVariable

##    def setMesh(self, newMesh):
##        newValues = self.getValue(points = newMesh.getCellCenters())
##        self.mesh = newMesh
##        self.setValue(newValues)
        
    def copy(self):
        
        return self.__class__(
            mesh = self.mesh, 
            name = self.name + "_old", 
            value = self.getValue(),
            hasOld = 0)
            
    def __call__(self, point = None):
        if point != None:
            return self[self.getMesh()._getNearestCellID(point)]
        else:
            return Variable.__call__(self)
##      return (self[i[0]] * self[i[1]] * (d[i[0]] + d[i[1]])) / (self[i[0]] * d[i[0]] + self[i[1]] * d[i[1]])
        
    def getValue(self, points = (), cells = ()):
        if points == () and cells == ():
            return Variable.getValue(self)
        elif cells != ():
            return numerix.take(Variable.getValue(self), [cell.getID() for cell in cells])
        else:
            return [self(point) for point in points]
        
    def setValue(self, value, cells = (), unit = None, where = None):
        """
        Patched values can be set by using either `cells` or `where`.

            >>> from fipy.meshes.grid1D import Grid1D
            >>> mesh = Grid1D(nx = 4)
            >>> v1 = CellVariable(value=(4,7,2,6), mesh=mesh)
            >>> print v1
            [ 4.  7.  2.  6.]
            >>> v1.setValue(4, where=(0, 0, 1, 1))
            >>> print v1
            [ 4.  7.  4.  4.]
            >>> v1.setValue((5,2,7,8), where=(0, 1, 1, 1))
            >>> print v1
            [ 4.  2.  7.  8.]
            >>> v1.setValue(3, unit = 'm')
            >>> print v1
            [ 3.  3.  3.  3.] m
            >>> import warnings
            >>> warnings.filterwarnings("ignore", "'where' should be used instead of 'cells'", DeprecationWarning)
            >>> v1.setValue(1, cells=mesh.getCells()[2:])
            >>> warnings.resetwarnings()
            >>> print v1
            [ 3.  3.  1.  1.] m
            >>> v1.setValue(4)
            >>> print v1
            [ 4.  4.  4.  4.]
            

        """
            
        if cells == ():
            Variable.setValue(self, value, unit = unit, where = where)
        else:
            import warnings
            warnings.warn("'where' should be used instead of 'cells'", DeprecationWarning, stacklevel=2)
            for cell in cells:
                self[cell.getID()] = value

            
    def getCellVolumeAverage(self):
        r"""
        Return the cell-volume-weighted average of the `CellVariable`:
            
        .. raw:: latex
        
           \[ <\phi>_\text{vol} 
           = \frac{\sum_\text{cells} \phi_\text{cell} V_\text{cell}}
               {\sum_\text{cells} V_\text{cell}} \]
        
        ..

            >>> from fipy.meshes.grid2D import Grid2D
            >>> mesh = Grid2D(nx = 3, ny = 1, dx = .5, dy = .1)
            >>> var = CellVariable(value = (1, 2, 6), mesh = mesh)
            >>> print var.getCellVolumeAverage()
            3.0
        """

        if self.volumeAverage is None:
            from cellVolumeAverageVariable import _CellVolumeAverageVariable
            self.volumeAverage = _CellVolumeAverageVariable(self)
        
        return self.volumeAverage

    def getGrad(self):
        r"""
        Return
        
        .. raw:: latex
        
           \( \nabla \phi \)
           
        as a `VectorCellVariable` (first-order gradient).
        """
        return self.getGaussGrad()

    def getGaussGrad(self):
        r"""
        Return

        .. raw:: latex

            $ \frac{1}{V_P} \sum_f \vec{n} \phi_f A_f $

        as a `VectorCellVariable` (first-order gradient).
            
        """
        if not hasattr(self, 'gaussGrad'):
            from gaussCellGradVariable import _GaussCellGradVariable
            self.gaussGrad = _GaussCellGradVariable(var = self, name = "%s_gauss_grad" % self.getName())
        
        return self.gaussGrad

    def getLeastSquaresGrad(self):
        r"""
        Return

        .. raw:: latex

            $\nabla \phi$, which is determined by solving for $\nabla
            \phi$ in the following matrix equation, $$ \nabla \phi \cdot \sum_f
            d_{AP}^2 \vec{n}_{AP} \otimes \vec{n}_{AP} = \sum_f
            d_{AP}^2 \left( \vec{n} \cdot \nabla \phi \right)_{AP} $$
            The matrix equation is derived by minimizing the following
            least squares sum, $$ F \left( \phi_x, \phi_y \right) =
            \sqrt{\sum_f \left( d_{AP} \vec{n}_{AP} \cdot \nabla \phi
            - d_{AP} \left( \vec{n}_{AP} \cdot \nabla \phi
            \right)_{AP} \right)^2 } $$

        Tests

            >>> from fipy import Grid2D
            >>> print CellVariable(mesh=Grid2D(nx=2, ny=2, dx=0.1, dy=2.0), value=(0,1,3,6)).getLeastSquaresGrad()
            [[8.0 1.2]
             [8.0 2.0]
             [24.0 1.2]
             [24.0 2.0]]

            >>> from fipy import Grid1D
            >>> print CellVariable(mesh=Grid1D(dx=(2.0, 1.0, 0.5)), value=(0, 1, 2)).getLeastSquaresGrad()
            [[0.461538461538]
             [0.8]
             [1.2]]

        """

        if not hasattr(self, 'leastSquaresGrad'):
            from leastSquaresCellGradVariable import _LeastSquaresCellGradVariable
            self.leastSquaresGrad = _LeastSquaresCellGradVariable(var = self, name = "%s_least_squares_grad" % self.getName())
        
        return self.leastSquaresGrad


    def getArithmeticFaceValue(self):
        r"""
        Returns a `FaceVariable` whose value corresponds to the arithmetic interpolation
        of the adjacent cells:
            
        .. raw:: latex
        
           \[ \phi_f = (\phi_1 - \phi_2) \frac{d_{f2}}{d_{12}} + \phi_2 \]
           
        ..
        
            >>> from fipy.meshes.grid1D import Grid1D
            >>> mesh = Grid1D(dx = (1., 1.))
            >>> var = CellVariable(mesh = mesh, value = (1, 2))
            >>> faceValue = var.getArithmeticFaceValue()[mesh.getInteriorFaces()[0]]
            >>> answer = (var[0] - var[1]) * (0.5 / 1.) + var[1]
            >>> numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
            1
            
            >>> mesh = Grid1D(dx = (2., 4.))
            >>> var = CellVariable(mesh = mesh, value = (1, 2))
            >>> faceValue = var.getArithmeticFaceValue()[mesh.getInteriorFaces()[0]]
            >>> answer = (var[0] - var[1]) * (1.0 / 3.0) + var[1]
            >>> numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
            1

            >>> mesh = Grid1D(dx = (10., 100.))
            >>> var = CellVariable(mesh = mesh, value = (1, 2))
            >>> faceValue = var.getArithmeticFaceValue()[mesh.getInteriorFaces()[0]]
            >>> answer = (var[0] - var[1]) * (5.0 / 55.0) + var[1]
            >>> numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
            1
        """
        if self.arithmeticFaceValue is None:
            from arithmeticCellToFaceVariable import _ArithmeticCellToFaceVariable
            self.arithmeticFaceValue = _ArithmeticCellToFaceVariable(self)

        return self.arithmeticFaceValue

    def getHarmonicFaceValue(self):
        r"""
        Returns a `FaceVariable` whose value corresponds to the harmonic interpolation
        of the adjacent cells:
            
        .. raw:: latex
        
           \[ \phi_f = \frac{\phi_1 \phi_2}{(\phi_2 - \phi_1) \frac{d_{f2}}{d_{12}} + \phi_1} \]
           
        ..
        
            >>> from fipy.meshes.grid1D import Grid1D
            >>> mesh = Grid1D(dx = (1., 1.))
            >>> var = CellVariable(mesh = mesh, value = (1, 2))
            >>> # faceValue = var.getHarmonicFaceValue()[mesh.getInteriorFaces()[0]]
            >>> faceValue = var.getHarmonicFaceValue()[mesh.getInteriorFaces()]
            >>> answer = var[0] * var[1] / ((var[1] - var[0]) * (0.5 / 1.) + var[0])
            >>> numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
            1
            
            >>> mesh = Grid1D(dx = (2., 4.))
            >>> var = CellVariable(mesh = mesh, value = (1, 2))
            >>> faceValue = var.getHarmonicFaceValue()[mesh.getInteriorFaces()[0]]
            >>> answer = var[0] * var[1] / ((var[1] - var[0]) * (1.0 / 3.0) + var[0])
            >>> numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
            1

            >>> mesh = Grid1D(dx = (10., 100.))
            >>> var = CellVariable(mesh = mesh, value = (1, 2))
            >>> faceValue = var.getHarmonicFaceValue()[mesh.getInteriorFaces()[0]]
            >>> answer = var[0] * var[1] / ((var[1] - var[0]) * (5.0 / 55.0) + var[0])
            >>> numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
            1
        """
        if self.harmonicFaceValue is None:
            from harmonicCellToFaceVariable import _HarmonicCellToFaceVariable
            self.harmonicFaceValue = _HarmonicCellToFaceVariable(self)

        return self.harmonicFaceValue

    def getFaceValue(self, func=None, **args):
        r"""
        Returns a `FaceVariable` whose value is evaluated by the function provided.

            >>> from fipy.meshes.grid1D import Grid1D
            >>> mesh = Grid1D(nx=2)
            >>> var = CellVariable(value=(0, 1), mesh=mesh)
            >>> print var.getFaceValue()
            [ 0.   0.5  1. ]
            >>> print var.getFaceValue(lambda a, b: numerix.minimum(a,b))
            [ 0.  0.  1.]
            
        """
        if func is None:
            return self.getArithmeticFaceValue()
        else:
            from callBackCellToFaceVariable import _CallBackCellToFaceVariable
            return _CallBackCellToFaceVariable(self, func, **args)

    def getFaceGrad(self):
        r"""
        Return
        
        .. raw:: latex
        
           \( \nabla \phi \)
           
        as a `VectorFaceVariable` using differencing for the normal direction(second-order gradient).
        """
        if self.faceGrad is None:
            from faceGradVariable import _FaceGradVariable
            self.faceGrad = _FaceGradVariable(self)

        return self.faceGrad

    def getFaceGradAverage(self):
        r"""
        Return
        
        .. raw:: latex
        
           \( \nabla \phi \)
           
        as a `VectorFaceVariable` using averaging for the normal direction(second-order gradient)
        """
        return self.getGrad().getArithmeticFaceValue()

    def getOld(self):
        """
        Return the values of the `CellVariable` from the previous
        solution sweep.

        Combinations of `CellVariable's` should also return old
        values.

            >>> from fipy.meshes.grid1D import Grid1D
            >>> mesh = Grid1D(nx = 2)
            >>> from fipy.variables.cellVariable import CellVariable
            >>> var1 = CellVariable(mesh = mesh, value = (2, 3), hasOld = 1)
            >>> var2 = CellVariable(mesh = mesh, value = (3, 4))
            >>> v = var1 * var2
            >>> print v
            [  6.  12.]
            >>> var1.setValue((3,2))
            >>> print v
            [ 9.  8.]
            >>> print v.getOld()
            [  6.  12.]

        The following small test is to correct for a bug when the
        operator does not just use variables.

            >>> v1 = var1 * 3
            >>> print v1
            [ 9.  6.]
            >>> print v1.getOld()
            [ 6.  9.]
            
        """
        if self.old is None:
            return self
        else:
            return self.old
##             import weakref
##          return weakref.proxy(self.old)

    def updateOld(self):
        """
        Set the values of the previous solution sweep to the current values.
        """
        if self.old is not None:
            self.old.setValue(self.getValue().copy())

    def _resetToOld(self):
        if self.old is not None:
            self.setValue(self.old.getValue())
            
    def _remesh(self, mesh):
        self.value = numerix.array(self.getValue(points = mesh.getCellCenters()))
        if self.old is not None:
            self.old._remesh(mesh)
        self.mesh = mesh
        self._markFresh()

    def _getShapeFromMesh(mesh):
        """
        Return the shape of this variable type, given a particular mesh.
        """
        return (mesh.getNumberOfCells(),)
    _getShapeFromMesh = staticmethod(_getShapeFromMesh)

    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base
        class for an operation result.
        """
        if other is None:
            return CellVariable
            
        # A CellVariable operating with a vector will produce a VectorCellVariable.
        # As a special case, if the number of cells equals the number of spatial dimensions,
        # treat tuple of that length as a vector, rather than as a list of scalars.
        if not isinstance(other, CellVariable) \
        and numerix.getShape(other) == (self.getMesh().getDim(),) \
        and (self.getMesh().getNumberOfCells(),) != (self.getMesh().getDim(),):
            from fipy.variables.vectorCellVariable import VectorCellVariable
            return VectorCellVariable
        else:
            return Variable._getArithmeticBaseClass(self, other)

    def _verifyShape(self, op, var0, var1, var0Array, var1Array, opShape, otherClass, rotateShape = True):
        from fipy.variables.vectorCellVariable import VectorCellVariable
        if isinstance(var1, VectorCellVariable) and self.getMesh() == var1.getMesh():
            return self._rotateShape(op, var1, var0, var1Array, var0Array, opShape)
        else:
            return Variable._verifyShape(self, op, var0, var1, var0Array, var1Array, opShape, otherClass, rotateShape)

##pickling
            
    def __getstate__(self):
        """
        Used internally to collect the necessary information to ``pickle`` the 
        `CellVariable` to persistent storage.
        """
        return {
            'mesh' : self.mesh,
            'name' : self.name,
            'value' : self.getValue(),
            'unit' : self.getUnit(),
            'old' : self.old
        }

    def __setstate__(self, dict):
        """
        Used internally to create a new `CellVariable` from ``pickled`` 
        persistent storage.
        """
        
        import sys
        self._refcount = sys.getrefcount(self)

        hasOld = 0
        if dict['old'] is not None:
            hasOld = 1

        self.__init__(mesh=dict['mesh'], name=dict['name'], value=dict['value'], unit=dict['unit'], hasOld=hasOld)
##         self.__init__(hasOld=hasOld, **dict)
        if self.old is not None:
            self.old.setValue(dict['old'].getValue())


class _ReMeshedCellVariable(CellVariable):
    def __init__(self, oldVar, newMesh):
        newValues = oldVar.getValue(points = newMesh.getCellCenters())
        CellVariable.__init__(self, newMesh, name = oldVar.name, value = newValues, unit = oldVar.getUnit())

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 

