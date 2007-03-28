#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "convectionTerm.py"
 #                                    created: 11/13/03 {11:39:03 AM} 
 #                                last update: 3/28/07 {10:17:10 AM} 
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-13 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.terms.faceTerm import FaceTerm
from fipy.variables.vectorFaceVariable import VectorFaceVariable
from fipy.variables.vectorCellVariable import VectorCellVariable
from fipy.variables.faceVariable import FaceVariable
from fipy.variables.cellVariable import CellVariable

from fipy.tools import numerix

class ConvectionTerm(FaceTerm):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, coeff = 1.0, diffusionTerm = None):
        """
        Create a `ConvectionTerm` object.
        
            >>> from fipy.meshes.grid1D import Grid1D
            >>> from fipy.variables.cellVariable import CellVariable
            >>> from fipy.variables.faceVariable import FaceVariable
            >>> from fipy.variables.vectorCellVariable import VectorCellVariable
            >>> from fipy.variables.vectorFaceVariable import VectorFaceVariable
            >>> m = Grid1D(nx = 2)
            >>> cv = CellVariable(mesh = m)
            >>> fv = FaceVariable(mesh = m)
            >>> vcv = VectorCellVariable(mesh = m)
            >>> vfv = VectorFaceVariable(mesh = m)
            >>> ConvectionTerm(coeff = cv)
            Traceback (most recent call last):
                ...
            TypeError: The coefficient must be a VectorFaceVariable, VectorCellVariable, or a vector value.
            >>> ConvectionTerm(coeff = fv)
            Traceback (most recent call last):
                ...
            TypeError: The coefficient must be a VectorFaceVariable, VectorCellVariable, or a vector value.
            >>> ConvectionTerm(coeff = vcv)
            ConvectionTerm(coeff = [[ 0.]
             [ 0.]
             [ 0.]])
            >>> ConvectionTerm(coeff = vfv)
            ConvectionTerm(coeff = [[ 0.]
             [ 0.]
             [ 0.]])
            >>> ConvectionTerm(coeff = (1,))
            ConvectionTerm(coeff=(1,))
            >>> from fipy.terms.explicitUpwindConvectionTerm import ExplicitUpwindConvectionTerm
            >>> ExplicitUpwindConvectionTerm(coeff = (0,)).solve(var = cv)
            >>> ExplicitUpwindConvectionTerm(coeff = 1).solve(var = cv)
            Traceback (most recent call last):
                ...
            TypeError: The coefficient must be a VectorFaceVariable, VectorCellVariable, or a vector value.

        
        :Parameters:
          - `coeff` : The `Term`'s coefficient value.
          - `diffusionTerm` : If a `DiffusionTerm` is given, the `ConvectionTerm` uses the diffusion coefficient to calculate the Peclet number.
        """
        self.diffusionTerm = diffusionTerm
        self.stencil = None
        
        if isinstance(coeff, VectorCellVariable):
            coeff = coeff.getArithmeticFaceValue()
            
        if isinstance(coeff, CellVariable) or isinstance(coeff, FaceVariable):
            raise TypeError, "The coefficient must be a VectorFaceVariable, VectorCellVariable, or a vector value."

        FaceTerm.__init__(self, coeff = coeff)
        
    def __neg__(self):
        """
        Negate the term.

           >>> -ConvectionTerm(coeff=1.0)
           ConvectionTerm(coeff=-1.0)
        """
        try:
            coeff = -self.coeff
        except:
            coeff = -numerix.array(self.coeff)

        return self.__class__(coeff=coeff, diffusionTerm=self.diffusionTerm)

    def _calcGeomCoeff(self, mesh):
        if not isinstance(self.coeff, VectorFaceVariable):
            self.coeff = VectorFaceVariable(mesh=mesh, value=self.coeff)
        
        projectedCoefficients = self.coeff * mesh._getOrientedAreaProjections()
        
        return projectedCoefficients.sum(1)
        
    def _getWeight(self, mesh):

        if self.stencil is None:

            if self.diffusionTerm is None:
                diffCoeff = 1e-20
            else:
                diffCoeff = self.diffusionTerm._getGeomCoeff(mesh)
                diffCoeff = diffCoeff * (diffCoeff != 0.) + 1e-20 * (diffCoeff == 0.)
                
##             alpha = self._Alpha(-self._getGeomCoeff(mesh) / diffCoeff)
            alpha = self._Alpha(self._diagonalSign * self._getGeomCoeff(mesh) / diffCoeff)
            
            self.stencil = {'implicit' : {'cell 1 diag'    : alpha,
                                          'cell 1 offdiag' : (1-alpha),
                                          'cell 2 diag'    : -(1-alpha),
                                          'cell 2 offdiag' : -alpha}}

        return self.stencil

    def _getDefaultSolver(self, solver):        
        if solver and not solver._canSolveAssymetric():
            import warnings
            warnings.warn("%s cannot solve assymetric matrices" % solver)
        from fipy.solvers.linearLUSolver import LinearLUSolver
        return solver or LinearLUSolver()

    def _verifyCoeffType(self, var):
        if not isinstance(self.coeff, VectorFaceVariable) \
        and numerix.getShape(self.coeff) != (var.getMesh().getDim(),):
            raise TypeError, "The coefficient must be a VectorFaceVariable, VectorCellVariable, or a vector value."

    def __add__(self, other):
        if isinstance(other, ConvectionTerm):
            if other.__class__ != self.__class__:
                raise TypeError, "ConvectionTerms must use the same scheme: %s != %s" % (self.__class__.__name__, other.__class__.__name__)
            return self.__class__(coeff=self.coeff + other.coeff)
        else:
            return FaceTerm.__add__(self, other)

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
