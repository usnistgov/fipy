#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "cellTerm.py"
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

from fipy.terms.nonDiffusionTerm import _NonDiffusionTerm
from fipy.tools import inline
from fipy.tools import numerix
from fipy.terms import AbstractBaseClassError
from fipy.variables.cellVariable import CellVariable
from fipy.variables.faceVariable import FaceVariable

__all__ = ["CellTerm"]

class CellTerm(_NonDiffusionTerm):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, coeff=1., var=None):
        if self.__class__ is CellTerm:
            raise AbstractBaseClassError

        from fipy.variables.variable import Variable

        if not isinstance(coeff, Variable):
            from fipy.variables.constant import _Constant
            coeff = _Constant(value=coeff)

        if isinstance(coeff, FaceVariable):
             raise TypeError, "The coefficient can not be a FaceVariable."

        _NonDiffusionTerm.__init__(self, coeff=coeff, var=var)
        self.coeffVectors = None
        self._var = None

    def _checkCoeff(self, var):
        if isinstance(self.coeff, CellVariable):
            shape = self.coeff.shape[:-1]
        else:
            shape = self.coeff.shape

        if var.rank == 1 and (shape == () or len(shape) == 1):
            if len(self.coeff.shape) == 2 and isinstance(self.coeff, CellVariable):
                self.coeff *= numerix.identity(var.shape[0])[...,numerix.newaxis]
            else:
                self.coeff *= numerix.identity(var.shape[0])
            if isinstance(self.coeff, CellVariable):
                shape = self.coeff.shape[:-1]
            else:
                shape = self.coeff.shape

        if var.rank == 0:
            if shape != ():
                raise TypeError, "The coefficient must be rank 0 for a rank 0 solution variable."

        if shape != () and len(shape) != 2 and shape[0] != shape[1]:
            raise TypeError, "The coefficient must be a rank-0 or rank-2 vector or a scalar value."

        if var.rank == 1:
            if shape == ():
                pass
            elif len(shape) != 2:
                raise TypeError, "The coefficient must be rank 2 or rank 0 for a rank 1 solution variable."
            elif var.shape[0] != shape[0]:
                raise TypeError, "The coefficient (N , N) shape must match the the solution variable (N,) shape."

    def _calcCoeffVectors_(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):
        coeff = self._getGeomCoeff(var)
        weight = self._getWeight(var, transientGeomCoeff, diffusionGeomCoeff)
        if hasattr(coeff, 'old'):
            old = coeff.old
        else:
            old = coeff

        self.coeffVectors = {
            'diagonal': coeff * weight['diagonal'],
            'old value': old * weight['old value'],
            'b vector': coeff * weight['b vector'],
            'new value': coeff * weight['new value']
        }

    def _getCoeffVectors_(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):
        if self.coeffVectors is None or var is not self._var:
            self._var = var
            self._calcCoeffVectors_(var=var, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

        return self.coeffVectors

    def _buildMatrixInline_(self, L, oldArray, b, dt, coeffVectors):
        oldArray = oldArray.value.ravel()
        N = len(oldArray)
        updatePyArray = numerix.zeros((N),'d')

        dt = self._checkDt(dt)

        inline._runInline("""
            b[i] += oldArray[i] * oldCoeff[i] / dt;
            b[i] += bCoeff[i];
            updatePyArray[i] += newCoeff[i] / dt;
            updatePyArray[i] += diagCoeff[i];
        """,b=b,
            oldArray=oldArray,
            oldCoeff=coeffVectors['old value'].ravel(),
            bCoeff=coeffVectors['b vector'].ravel(),
            newCoeff=coeffVectors['new value'].ravel(),
            diagCoeff=coeffVectors['diagonal'].ravel(),
            updatePyArray=updatePyArray,
            ni=len(updatePyArray),
            dt=dt)

        L.addAtDiagonal(updatePyArray)

    def _buildMatrixNoInline_(self, L, oldArray, b, dt, coeffVectors):
        ids = self._reshapeIDs(oldArray, numerix.arange(oldArray.shape[-1]))
        b += (oldArray.value[numerix.newaxis] * coeffVectors['old value']).sum(-2).ravel() / dt
        b += coeffVectors['b vector'][numerix.newaxis].sum(-2).ravel()
        L.addAt(coeffVectors['new value'].ravel() / dt, ids.ravel(), ids.swapaxes(0,1).ravel())
        L.addAt(coeffVectors['diagonal'].ravel(), ids.ravel(), ids.swapaxes(0,1).ravel())

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):

        b = numerix.zeros(var.shape,'d').ravel()
        L = SparseMatrix(mesh=var.mesh)

        coeffVectors = self._getCoeffVectors_(var=var, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

        dt = self._checkDt(dt)

        if inline.doInline and var.rank == 0:
            self._buildMatrixInline_(L=L, oldArray=var.old, b=b, dt=dt, coeffVectors=coeffVectors)
        else:
            self._buildMatrixNoInline_(L=L, oldArray=var.old, b=b, dt=dt, coeffVectors=coeffVectors)

        return (var, L, b)

    def _test(self):
        """
        The following tests demonstrate how the `CellVariable` objects
        interact with other types of `Variable` objects.

            >>> from fipy.meshes import Grid1D
            >>> from fipy.variables.cellVariable import CellVariable
            >>> from fipy.variables.faceVariable import FaceVariable
            >>> from fipy.terms.transientTerm import TransientTerm
            >>> m = Grid1D(nx=2)
            >>> cv = CellVariable(mesh=m)
            >>> fv = FaceVariable(mesh=m)
            >>> vcv = CellVariable(mesh=m, rank=1)
            >>> vfv = FaceVariable(mesh=m, rank=1)

            >>> __CellTerm(coeff=cv)
            __CellTerm(coeff=CellVariable(value=array([ 0.,  0.]), mesh=UniformGrid1D(dx=1.0, nx=2)))
            >>> __CellTerm(coeff=1)
            __CellTerm(coeff=1)
            >>> __CellTerm(coeff=fv)
            Traceback (most recent call last):
                ...
            TypeError: The coefficient can not be a FaceVariable.
            >>> TransientTerm(coeff=vcv).solve(cv)
            Traceback (most recent call last):
                ...
            TypeError: The coefficient must be rank 0 for a rank 0 solution variable.
            >>> __CellTerm(coeff=vfv)
            Traceback (most recent call last):
                ...
            TypeError: The coefficient can not be a FaceVariable.
            >>> TransientTerm(coeff=(1,)).solve(cv)
            Traceback (most recent call last):
                ...
            TypeError: The coefficient must be rank 0 for a rank 0 solution variable.

        """
        pass

class __CellTerm(CellTerm):
    """
    Dummy subclass for tests
    """
    pass


def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
