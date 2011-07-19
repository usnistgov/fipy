#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "baseConvectionTerm.py"
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
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.terms.faceTerm import FaceTerm
from fipy.variables.meshVariable import _MeshVariable
from fipy.variables.faceVariable import FaceVariable
from fipy.variables.cellVariable import CellVariable
from fipy.terms import AbstractBaseClassError
from fipy.terms import VectorCoeffError

from fipy.tools import numerix

class _BaseConvectionTerm(FaceTerm):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, coeff=1.0, diffusionTerm=None, var=None):
        """
        Create a `_BaseConvectionTerm` object.
        
            >>> from fipy import *
            >>> m = Grid1D(nx = 2)
            >>> cv = CellVariable(mesh = m)
            >>> fv = FaceVariable(mesh = m)
            >>> vcv = CellVariable(mesh=m, rank=1)
            >>> vfv = FaceVariable(mesh=m, rank=1)
            >>> __ConvectionTerm(coeff = cv)
            Traceback (most recent call last):
                ...
            VectorCoeffError: The coefficient must be a vector value.
            >>> __ConvectionTerm(coeff = fv)
            Traceback (most recent call last):
                ...
            VectorCoeffError: The coefficient must be a vector value.
            >>> __ConvectionTerm(coeff = vcv)
            __ConvectionTerm(coeff=_ArithmeticCellToFaceVariable(value=array([[ 0.,  0.,  0.]]), mesh=UniformGrid1D(dx=1.0, nx=2)))
            >>> __ConvectionTerm(coeff = vfv)
            __ConvectionTerm(coeff=FaceVariable(value=array([[ 0.,  0.,  0.]]), mesh=UniformGrid1D(dx=1.0, nx=2)))
            >>> __ConvectionTerm(coeff = (1,))
            __ConvectionTerm(coeff=(1,))
            >>> ExplicitUpwindConvectionTerm(coeff = (0,)).solve(var=cv, solver=DummySolver())
            >>> ExplicitUpwindConvectionTerm(coeff = 1).solve(var=cv, solver=DummySolver())
            Traceback (most recent call last):
                ...
            VectorCoeffError: The coefficient must be a vector value.
            >>> m2 = Grid2D(nx=2, ny=1)
            >>> cv2 = CellVariable(mesh=m2)
            >>> vcv2 = CellVariable(mesh=m2, rank=1)
            >>> vfv2 = FaceVariable(mesh=m2, rank=1)
            >>> __ConvectionTerm(coeff=vcv2)
            __ConvectionTerm(coeff=_ArithmeticCellToFaceVariable(value=array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.],
                   [ 0.,  0.,  0.,  0.,  0.,  0.,  0.]]), mesh=UniformGrid2D(dx=1.0, dy=1.0, nx=2, ny=1)))
            >>> __ConvectionTerm(coeff=vfv2)
            __ConvectionTerm(coeff=FaceVariable(value=array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.],
                   [ 0.,  0.,  0.,  0.,  0.,  0.,  0.]]), mesh=UniformGrid2D(dx=1.0, dy=1.0, nx=2, ny=1)))
            >>> ExplicitUpwindConvectionTerm(coeff = ((0,),(0,))).solve(var=cv2, solver=DummySolver())
            >>> ExplicitUpwindConvectionTerm(coeff = (0,0)).solve(var=cv2, solver=DummySolver())

        
        :Parameters:
          - `coeff` : The `Term`'s coefficient value.
          - `diffusionTerm` : **deprecated**. The Peclet number is calculated automatically.
        """
        if self.__class__ is _BaseConvectionTerm:
            raise AbstractBaseClassError
            
        if diffusionTerm is not None:
            import warnings
            warnings.warn("The Peclet number is calculated automatically. diffusionTerm will be ignored.", DeprecationWarning, stacklevel=2)

        self.stencil = None
        
        if isinstance(coeff, _MeshVariable) and coeff.rank < 1:
            raise VectorCoeffError

        if isinstance(coeff, CellVariable):
            coeff = coeff.arithmeticFaceValue

        FaceTerm.__init__(self, coeff=coeff, var=var)
        
    def _calcGeomCoeff(self, var):
        mesh = var.mesh

        if not isinstance(self.coeff, FaceVariable):
            shape = numerix.array(self.coeff).shape
            if shape != () and shape[-1] == 1:
                shape = shape[:-1]
            
            self.coeff = FaceVariable(mesh=mesh, elementshape=shape, value=self.coeff)

        projectedCoefficients = self.coeff * mesh._orientedAreaProjections

        return projectedCoefficients.sum(0)
        
    def _getWeight(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):

        if self.stencil is None:

            small = -1e-20

            if diffusionGeomCoeff is None:
                diffCoeff = small
            else:
                diffCoeff = diffusionGeomCoeff[0]
                if diffCoeff is None:
                    diffCoeff = small
                else:
                    diffCoeff = diffCoeff.numericValue
                    diffCoeff = (diffCoeff == 0) * small + diffCoeff

            alpha = self._alpha(-self._getGeomCoeff(var) / diffCoeff)
            
            self.stencil = {'implicit' : {'cell 1 diag'    : alpha,
                                          'cell 1 offdiag' : (1-alpha),
                                          'cell 2 diag'    : -(1-alpha),
                                          'cell 2 offdiag' : -alpha}}

        return self.stencil

    def _checkVar(self, var):
        FaceTerm._checkVar(self, var)
        
        if not (isinstance(self.coeff, FaceVariable) and self.coeff.rank == 1) \
        and numerix.getShape(self.coeff) != (var.mesh.dim,):
            raise VectorCoeffError

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=1., transientGeomCoeff=None, diffusionGeomCoeff=None):
        
        var, L, b = FaceTerm._buildMatrix(self, var, SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

        mesh = var.mesh

        if (not hasattr(self, 'constraintL')) or (not hasattr(self, 'constraintB')):

            constraintMask = var.faceGrad.constraintMask | var.arithmeticFaceValue.constraintMask

            weight = self._getWeight(var, transientGeomCoeff, diffusionGeomCoeff)

            if weight.has_key('implicit'):
                alpha = weight['implicit']['cell 1 diag']
            else:
                alpha = 0.0

            exteriorCoeff =  self.coeff.faceValue * mesh.exteriorFaces

            self.constraintL = (alpha * constraintMask * exteriorCoeff).divergence * mesh.cellVolumes
            self.constraintB =  -((1 - alpha) * var.arithmeticFaceValue * constraintMask * exteriorCoeff).divergence * mesh.cellVolumes

        ids = self._reshapeIDs(var, numerix.arange(mesh.numberOfCells))
        L.addAt(numerix.array(self.constraintL).ravel(), ids.ravel(), ids.swapaxes(0,1).ravel())
        b += numerix.reshape(self.constraintB.value, ids.shape).sum(0).ravel()

        return (var, L, b)

class __ConvectionTerm(_BaseConvectionTerm): 
    """
    Dummy subclass for tests
    """
    pass 

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
