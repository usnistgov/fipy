#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "convectionTerm.py"
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
from fipy.solvers import DefaultAsymmetricSolver

from fipy.tools import numerix

class ConvectionTerm(FaceTerm):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, coeff=1.0, diffusionTerm=None):
        """
        Create a `ConvectionTerm` object.
        
            >>> from fipy.meshes.grid1D import Grid1D
            >>> from fipy.variables.cellVariable import CellVariable
            >>> from fipy.variables.faceVariable import FaceVariable
            >>> m = Grid1D(nx = 2)
            >>> cv = CellVariable(mesh = m)
            >>> fv = FaceVariable(mesh = m)
            >>> vcv = CellVariable(mesh=m, rank=1)
            >>> vfv = FaceVariable(mesh=m, rank=1)
            >>> __ConvectionTerm(coeff = cv)
            Traceback (most recent call last):
                ...
            TypeError: The coefficient must be a vector value.
            >>> __ConvectionTerm(coeff = fv)
            Traceback (most recent call last):
                ...
            TypeError: The coefficient must be a vector value.
            >>> __ConvectionTerm(coeff = vcv)
            __ConvectionTerm(coeff=_ArithmeticCellToFaceVariable(value=array([[ 0.,  0.,  0.]]), mesh=UniformGrid1D(dx=1.0, nx=2)))
            >>> __ConvectionTerm(coeff = vfv)
            __ConvectionTerm(coeff=FaceVariable(value=array([[ 0.,  0.,  0.]]), mesh=UniformGrid1D(dx=1.0, nx=2)))
            >>> __ConvectionTerm(coeff = (1,))
            __ConvectionTerm(coeff=(1,))
            >>> from fipy.terms.explicitUpwindConvectionTerm import ExplicitUpwindConvectionTerm
            >>> ExplicitUpwindConvectionTerm(coeff = (0,)).solve(var = cv)
            >>> ExplicitUpwindConvectionTerm(coeff = 1).solve(var = cv)
            Traceback (most recent call last):
                ...
            TypeError: The coefficient must be a vector value.
            >>> from fipy.meshes.grid2D import Grid2D
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
            >>> ExplicitUpwindConvectionTerm(coeff = ((0,),(0,))).solve(var=cv2)
            >>> ExplicitUpwindConvectionTerm(coeff = (0,0)).solve(var=cv2)

        
        :Parameters:
          - `coeff` : The `Term`'s coefficient value.
          - `diffusionTerm` : **deprecated**. The Peclet number is calculated automatically.
        """
        if self.__class__ is ConvectionTerm:
            raise NotImplementedError, "can't instantiate abstract base class"
            
        if diffusionTerm is not None:
            import warnings
            warnings.warn("The Peclet number is calculated automatically. diffusionTerm will be ignored.", DeprecationWarning, stacklevel=2)

        self.stencil = None
        
        if isinstance(coeff, _MeshVariable) and coeff.getRank() != 1:
            raise TypeError, "The coefficient must be a vector value."

        if isinstance(coeff, CellVariable):
            coeff = coeff.getArithmeticFaceValue()

        FaceTerm.__init__(self, coeff = coeff)
        
    def _calcGeomCoeff(self, mesh):
        if not isinstance(self.coeff, FaceVariable):
            self.coeff = FaceVariable(mesh=mesh, value=self.coeff, rank=1)
        
        projectedCoefficients = self.coeff * mesh._getOrientedAreaProjections()
        
        return projectedCoefficients.sum(0)
        
    def _getWeight(self, mesh, equation=None):

        if self.stencil is None:

            small = -1e-20
            
            if equation is None:
                diffCoeff = small
            else:
                diffCoeff = equation._getDiffusiveGeomCoeff(mesh)
                if diffCoeff is None:
                    diffCoeff = small
                else:
                    diffCoeff = diffCoeff.getNumericValue()
                    diffCoeff = (diffCoeff == 0) * small + diffCoeff

            alpha = self._Alpha(-self._getGeomCoeff(mesh) / diffCoeff)
            
            self.stencil = {'implicit' : {'cell 1 diag'    : alpha,
                                          'cell 1 offdiag' : (1-alpha),
                                          'cell 2 diag'    : -(1-alpha),
                                          'cell 2 offdiag' : -alpha}}

        return self.stencil

    def _getDefaultSolver(self, solver, *args, **kwargs):        
        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve assymetric matrices" % solver)
        return solver or DefaultAsymmetricSolver(*args, **kwargs)

    def _verifyCoeffType(self, var):
        if not (isinstance(self.coeff, FaceVariable) and self.coeff.getRank() == 1) \
        and numerix.getShape(self.coeff) != (var.getMesh().getDim(),):
            raise TypeError, "The coefficient must be a vector value."

    def __add__(self, other):
        if isinstance(other, ConvectionTerm):
            if other.__class__ != self.__class__:
                raise TypeError, "ConvectionTerms must use the same scheme: %s != %s" % (self.__class__.__name__, other.__class__.__name__)
            return self.__class__(coeff=self.coeff + other.coeff)
        else:
            return FaceTerm.__add__(self, other)

class __ConvectionTerm(ConvectionTerm): 
    """
    Dummy subclass for tests
    """
    pass 

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
