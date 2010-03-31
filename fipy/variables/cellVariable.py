#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cellVariable.py"
 #
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

from fipy.variables.meshVariable import _MeshVariable
from fipy.tools import numerix

class CellVariable(_MeshVariable):
    """
    Represents the field of values of a variable on a `Mesh`.

    A `CellVariable` can be ``pickled`` to persistent storage (disk) for later use:
        
        >>> from fipy.meshes.grid2D import Grid2D
        >>> mesh = Grid2D(dx = 1., dy = 1., nx = 10, ny = 10)
        
        >>> var = CellVariable(mesh = mesh, value = 1., hasOld = 1, name = 'test')
        >>> x, y = mesh.getCellCenters()
        >>> var.setValue(x * y)

        >>> from fipy.tools import dump        
        >>> (f, filename) = dump.write(var, extension = '.gz')
        >>> unPickledVar = dump.read(filename, f)
        
        >>> print var.allclose(unPickledVar, atol = 1e-10, rtol = 1e-10)
        1
        
    """

    def __init__(self, mesh, name='', value=0., rank=None, elementshape=None, unit=None, hasOld=0):
        _MeshVariable.__init__(self, mesh=mesh, name=name, value=value, 
                               rank=rank, elementshape=elementshape, unit=unit)

        if hasOld:
            self.old = self.copy()
        else:
            self.old = None
            
    def _getVariableClass(self):
        return CellVariable
        
    def _OperatorVariableClass(self, baseClass=None):
        """
            >>> from fipy.meshes.grid1D import Grid1D

            >>> a = CellVariable(mesh=Grid1D(nx=1), value=1)

            >>> c = -a
            >>> b = c.getOld() + 3
            >>> print b
            [2]
            >>> b.getsctype() == numerix.NUMERIX.obj2sctype(numerix.array(1))
            True
            
        replacing with the same thing is no problem
        
            >>> a.setValue(3)
            >>> b = c.getOld() + 3
            >>> print b
            [0]
            
        replacing with multiple copies causes the reference counting problem
        
            >>> a.setValue(3)
            >>> b = (c + c).getOld() + 3
            >>> print b
            [-3]
            
        the order matters
        
            >>> b = (c + c).getOld() + 3
            >>> a.setValue(2)
            >>> print b
            [-1]
        """
        baseClass = _MeshVariable._OperatorVariableClass(self, 
                                                         baseClass=baseClass)
                                     
        class _CellOperatorVariable(baseClass):
            def getOld(self):
                if self.old is None:
                    oldVar = []
                    for v in self.var:
                        if hasattr(v, "getOld"):
                            oldVar.append(v.getOld())
                        else:
                            oldVar.append(v)
                    
                    self.old = self.__class__(op=self.op, var=oldVar, 
                                              opShape=self.opShape, 
                                              canInline=self.canInline)
                                  
                return self.old
                
        return _CellOperatorVariable
        
    def copy(self):
        
        return self._getArithmeticBaseClass()(mesh=self.mesh, 
                                              name=self.name + "_old", 
                                              value=self.getValue(),
                                              hasOld=False)
                
    def _getGlobalNumberOfElements(self):
        return self.mesh.globalNumberOfCells
        
    def _getGlobalOverlappingIDs(self):
        return self.mesh._getGlobalOverlappingCellIDs()

    def _getLocalNonOverlappingIDs(self):
        return self.mesh._getLocalNonOverlappingCellIDs()

    def getGlobalValue(self):
        """Concatenate and return values from all processors
        
        When running on a single processor, the result is identical to
        :meth:`~fipy.variables.variable.Variable.getValue`.
        """
        return self._getGlobalValue(self.mesh._getLocalNonOverlappingCellIDs(), 
                                    self.mesh._getGlobalNonOverlappingCellIDs())

    def setValue(self, value, unit = None, where = None):
        _MeshVariable.setValue(self, value=self._globalToLocalValue(value), unit=unit, where=where)

    def __call__(self, points=None, order=0, nearestCellIDs=None):
        r"""
        Interpolates the CellVariable to a set of points using a
        method that has a memory requirement on the order of Ncells by
        Npoints in general, but uses only Ncells when the
        CellVariable's mesh is a UniformGrid object.

        :Parameters:

           - `points`: A point or set of points in the format (X, Y, Z)
           - `order`: The order of interpolation, 0 or 1, default is 0 
           - `nearestCellIDs` : Optional argument if user can calculate own
             nearest cell IDs array, shape should be same as points

        Tests

            >>> from fipy import *
            >>> m = Grid2D(nx=3, ny=2)
            >>> v = CellVariable(mesh=m, value=m.getCellCenters()[0])
            >>> print v(((0., 1.1, 1.2), (0., 1., 1.)))
            [ 0.5  1.5  1.5]
            >>> print v(((0., 1.1, 1.2), (0., 1., 1.)), order=1)
            [ 0.25  1.1   1.2 ]
            >>> m0 = Grid2D(nx=2, ny=2, dx=1., dy=1.)
            >>> m1 = Grid2D(nx=4, ny=4, dx=.5, dy=.5)
            >>> x, y = m0.getCellCenters()
            >>> v0 = CellVariable(mesh=m0, value=x * y)
            >>> print v0(m1.getCellCenters().getGlobalValue())
            [ 0.25  0.25  0.75  0.75  0.25  0.25  0.75  0.75  0.75  0.75  2.25  2.25
              0.75  0.75  2.25  2.25]
            >>> print v0(m1.getCellCenters().getGlobalValue(), order=1)
            [ 0.125  0.25   0.5    0.625  0.25   0.375  0.875  1.     0.5    0.875
              1.875  2.25   0.625  1.     2.25   2.625]

        """           
        if points is not None:

            if nearestCellIDs is None:
                nearestCellIDs = self.getMesh()._getNearestCellID(points)

            if order == 0:
                return self.getGlobalValue()[..., nearestCellIDs]

            elif order == 1:
                ##cellID = self.getMesh()._getNearestCellID(points)
##                return self[...,self.getMesh()._getNearestCellID(points)] + numerix.dot(points - self.getMesh().getCellCenters()[...,cellID], self.getGrad()[...,cellID])
                return (self.getGlobalValue()[..., nearestCellIDs] 
                        + numerix.dot(points - self.getMesh().getCellCenters().getGlobalValue()[...,nearestCellIDs], 
                                      self.getGrad().getGlobalValue()[...,nearestCellIDs]))

            else:
                raise ValueError, 'order should be either 0 or 1'

        else:
            return _MeshVariable.__call__(self)
        
    def getCellVolumeAverage(self):
        r"""
        Return the cell-volume-weighted average of the `CellVariable`:
            
        .. math::
        
           <\phi>_\text{vol} 
           = \frac{\sum_\text{cells} \phi_\text{cell} V_\text{cell}}
               {\sum_\text{cells} V_\text{cell}}
        
        >>> from fipy.meshes.grid2D import Grid2D
        >>> from fipy.variables.cellVariable import CellVariable
        >>> mesh = Grid2D(nx = 3, ny = 1, dx = .5, dy = .1)
        >>> var = CellVariable(value = (1, 2, 6), mesh = mesh)
        >>> print var.getCellVolumeAverage()
        3.0
        """

        if not hasattr(self, 'volumeAverage'):
            from cellVolumeAverageVariable import _CellVolumeAverageVariable
            self.volumeAverage = _CellVolumeAverageVariable(self)
        
        return self.volumeAverage

    def getGrad(self):
        r"""
        Return :math:`\nabla \phi` as a rank-1 `CellVariable` (first-order
        gradient).
        """
        return self.getGaussGrad()

    def getGaussGrad(self):
        r"""
        Return :math:`\frac{1}{V_P} \sum_f \vec{n} \phi_f A_f`
        as a rank-1 `CellVariable` (first-order gradient).
            
        """
        if not hasattr(self, 'gaussGrad'):
            from gaussCellGradVariable import _GaussCellGradVariable
            self.gaussGrad = _GaussCellGradVariable(var = self, name = "%s_gauss_grad" % self.getName())
        
        return self.gaussGrad

    def getLeastSquaresGrad(self):
        r"""
        Return :math:`\nabla \phi`, which is determined by solving for :math:`\nabla \phi`
        in the following matrix equation, 
        
        .. math::
            
           \nabla \phi \cdot \sum_f d_{AP}^2 \vec{n}_{AP} \otimes \vec{n}_{AP} =
           \sum_f d_{AP}^2 \left( \vec{n} \cdot \nabla \phi \right)_{AP}
        
        The matrix equation is derived by minimizing
        the following least squares sum, 
        
        .. math:: 
           
           F \left( \phi_x, \phi_y \right) = \sqrt{\sum_f \left( d_{AP}
           \vec{n}_{AP} \cdot \nabla \phi - d_{AP} \left( \vec{n}_{AP} \cdot
           \nabla \phi \right)_{AP} \right)^2 }

        Tests

        >>> from fipy import Grid2D
        >>> m = Grid2D(nx=2, ny=2, dx=0.1, dy=2.0)
        >>> print numerix.allclose(CellVariable(mesh=m, value=(0,1,3,6)).getLeastSquaresGrad().getGlobalValue(), \
        ...                                     [[8.0, 8.0, 24.0, 24.0],
        ...                                      [1.2, 2.0, 1.2, 2.0]])
        True

        >>> from fipy import Grid1D
        >>> print numerix.allclose(CellVariable(mesh=Grid1D(dx=(2.0, 1.0, 0.5)), 
        ...                                     value=(0, 1, 2)).getLeastSquaresGrad().getGlobalValue(), [[0.461538461538, 0.8, 1.2]])
        True
        """

        if not hasattr(self, 'leastSquaresGrad'):
            from leastSquaresCellGradVariable import _LeastSquaresCellGradVariable
            self.leastSquaresGrad = _LeastSquaresCellGradVariable(var = self, name = "%s_least_squares_grad" % self.getName())
        
        return self.leastSquaresGrad


    def getArithmeticFaceValue(self):
        r"""
        Returns a `FaceVariable` whose value corresponds to the arithmetic interpolation
        of the adjacent cells:
            
        .. math::
        
           \phi_f = (\phi_1 - \phi_2) \frac{d_{f2}}{d_{12}} + \phi_2
           
        >>> from fipy.meshes.grid1D import Grid1D
        >>> from fipy import numerix
        >>> mesh = Grid1D(dx = (1., 1.))
        >>> L = 1
        >>> R = 2
        >>> var = CellVariable(mesh = mesh, value = (L, R))
        >>> faceValue = var.getArithmeticFaceValue()[mesh.getInteriorFaces().getValue()]
        >>> answer = (R - L) * (0.5 / 1.) + L
        >>> print numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
        True
        
        >>> mesh = Grid1D(dx = (2., 4.))
        >>> var = CellVariable(mesh = mesh, value = (L, R))
        >>> faceValue = var.getArithmeticFaceValue()[mesh.getInteriorFaces().getValue()]
        >>> answer = (R - L) * (1.0 / 3.0) + L
        >>> print numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
        True

        >>> mesh = Grid1D(dx = (10., 100.))
        >>> var = CellVariable(mesh = mesh, value = (L, R))
        >>> faceValue = var.getArithmeticFaceValue()[mesh.getInteriorFaces().getValue()]
        >>> answer = (R - L) * (5.0 / 55.0) + L
        >>> print numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
        True
        """
        if not hasattr(self, 'arithmeticFaceValue'):
            from arithmeticCellToFaceVariable import _ArithmeticCellToFaceVariable
            self.arithmeticFaceValue = _ArithmeticCellToFaceVariable(self)

        return self.arithmeticFaceValue

    getFaceValue = getArithmeticFaceValue

    def getMinmodFaceValue(self):
        r"""
        Returns a `FaceVariable` with a value that is the minimum of
        the absolute values of the adjacent cells. If the values are
        of opposite sign then the result is zero:
            
        .. math::
        
           \phi_f = \begin{cases}
                          \phi_1& \text{when $|\phi_1| \le |\phi_2|$},\\
                          \phi_2& \text{when $|\phi_2| < |\phi_1|$},\\
                          0 & \text{when $\phi1 \phi2 < 0$}
                    \end{cases}
                       
        >>> from fipy import *
        >>> print CellVariable(mesh=Grid1D(nx=2), value=(1, 2)).getMinmodFaceValue()
        [1 1 2]
        >>> print CellVariable(mesh=Grid1D(nx=2), value=(-1, -2)).getMinmodFaceValue()
        [-1 -1 -2]
        >>> print CellVariable(mesh=Grid1D(nx=2), value=(-1, 2)).getMinmodFaceValue()
        [-1  0  2]
        """
        if not hasattr(self, 'minmodFaceValue'):
            from minmodCellToFaceVariable import _MinmodCellToFaceVariable
            self.minmodFaceValue = _MinmodCellToFaceVariable(self)

        return self.minmodFaceValue

    def getHarmonicFaceValue(self):
        r"""
        Returns a `FaceVariable` whose value corresponds to the harmonic interpolation
        of the adjacent cells:
            
        .. math::
        
           \phi_f = \frac{\phi_1 \phi_2}{(\phi_2 - \phi_1) \frac{d_{f2}}{d_{12}} + \phi_1}
           
        >>> from fipy.meshes.grid1D import Grid1D
        >>> from fipy import numerix
        >>> mesh = Grid1D(dx = (1., 1.))
        >>> L = 1
        >>> R = 2
        >>> var = CellVariable(mesh = mesh, value = (L, R))
        >>> faceValue = var.getHarmonicFaceValue()[mesh.getInteriorFaces().getValue()]
        >>> answer = L * R / ((R - L) * (0.5 / 1.) + L)
        >>> print numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
        True
        
        >>> mesh = Grid1D(dx = (2., 4.))
        >>> var = CellVariable(mesh = mesh, value = (L, R))
        >>> faceValue = var.getHarmonicFaceValue()[mesh.getInteriorFaces().getValue()]
        >>> answer = L * R / ((R - L) * (1.0 / 3.0) + L)
        >>> print numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
        True

        >>> mesh = Grid1D(dx = (10., 100.))
        >>> var = CellVariable(mesh = mesh, value = (L, R))
        >>> faceValue = var.getHarmonicFaceValue()[mesh.getInteriorFaces().getValue()]
        >>> answer = L * R / ((R - L) * (5.0 / 55.0) + L)
        >>> print numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
        True
        """
        if not hasattr(self, 'harmonicFaceValue'):
            from harmonicCellToFaceVariable import _HarmonicCellToFaceVariable
            self.harmonicFaceValue = _HarmonicCellToFaceVariable(self)

        return self.harmonicFaceValue

    def getFaceGrad(self):
        r"""
        Return :math:`\nabla \phi` as a rank-1 `FaceVariable` using differencing
        for the normal direction(second-order gradient).
        """
        if not hasattr(self, 'faceGrad'):
            from faceGradVariable import _FaceGradVariable
            self.faceGrad = _FaceGradVariable(self)

        return self.faceGrad

    def getFaceGradAverage(self):
        r"""
        Return :math:`\nabla \phi` as a rank-1 `FaceVariable` using averaging
        for the normal direction(second-order gradient)
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
        [ 6 12]
        >>> var1.setValue((3,2))
        >>> print v
        [9 8]
        >>> print v.getOld()
        [ 6 12]

        The following small test is to correct for a bug when the
        operator does not just use variables.

        >>> v1 = var1 * 3
        >>> print v1
        [9 6]
        >>> print v1.getOld()
        [6 9]
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
            
        return _MeshVariable._getArithmeticBaseClass(self, other)

##pickling
            
    def __getstate__(self):
        """
        Used internally to collect the necessary information to ``pickle`` the 
        `CellVariable` to persistent storage.
        """
        return {
            'mesh' : self.mesh,
            'name' : self.name,
            'value' : self.getGlobalValue(),
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

