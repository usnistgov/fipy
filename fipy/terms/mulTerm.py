#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "diffusionTerm.py"
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

from fipy.terms.term import Term
from fipy.tools import numerix

class _MulTerm(Term):

    r"""
    This term represents a CellVariable mutiplied by a term or equation.

    """

    uniqueID = 0

    def __init__(self, coeff=1., term=None):
        """
        Create a `_MulTerm`.

        :Parameters:
          - `coeff`: Float` or `CellVariable`.
          - `term`: A subclass of `Term`
          
        """

        Term.__init__(self, coeff=coeff)
        self.term = term

        self.uniqueID += 1
        
    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=1., equation=None):

        L , b = self.term._buildMatrix(var,
                                       SparseMatrix,
                                       boundaryConditions=boundaryConditions,
                                       dt=dt,
                                       equation=equation)

        if len(numerix.shape(self.coeff)) == 0:
            matcoeff = float(self.coeff)
        else:
            matcoeff = SparseMatrix(mesh=var.getMesh(), bandwidth=1)
            matcoeff.addAtDiagonal(self.coeff)

        return (matcoeff * L, self.coeff * b)

    def _getDefaultSolver(self, solver, *args, **kwargs):

        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve assymetric matrices" % solver)

        from fipy.solvers import DefaultAsymmetricSolver
        return solver or DefaultAsymmetricSolver(*args, **kwargs)

    def __repr__(self):
        """
        The representation of a `_MulTerm` object is given by,

           >>> from fipy import *
           >>> print _MulTerm(coeff=1., term=DiffusionTerm())
           1.0 * DiffusionTerm(coeff=(1.0,))
           
        """
        return "%s * %s" % (repr(self.coeff), repr(self.term))

    def __neg__(self):
        return self.__class__(coeff=self.coeff, term=-self.term)

    def __add__(self, other):
        if self.__class__ == other.__class__:
            return self._add(other)
        else:
            return Term.__add__(self, other)

    def _isAdditive(self):
        return False

    def _test(self):
        r"""
        Test stuff

           >>> from fipy import *
           >>> L = 1.
           >>> nx = 100
           >>> m = Grid1D(nx=nx, dx=L / nx)
           >>> v = CellVariable(mesh=m, value=1.)
           >>> eqn = DiffusionTerm() / v - 1
           >>> BCs = (FixedValue(faces=m.getFacesLeft(), value=0.),
           ...        FixedValue(faces=m.getFacesRight(), value=1.))
           >>> res = 1.
           >>> sweep = 0
           >>> while res > 1e-8 and sweep < 100:
           ...     res = eqn.sweep(v, boundaryConditions=BCs)
           ...     sweep += 1
           >>> x = m.getCellCenters()[0]
           >>> answer = (numerix.exp(x) - numerix.exp(-x)) / (numerix.exp(L) - numerix.exp(-L))
           >>> print numerix.allclose(v, answer, rtol=2e-5)
           True

           >>> v.setValue(0.)
           >>> eqn = DiffusionTerm(0.2) * 5. - 5. * ImplicitSourceTerm(0.2)
           >>> eqn.solve(v, boundaryConditions=BCs)
           >>> print numerix.allclose(v, answer, rtol=2e-5)
           True

           >>> v.setValue(0.)
           >>> eqn = 2. * (DiffusionTerm(1.) - ImplicitSourceTerm(.5)) - DiffusionTerm(1.)
           >>> eqn.solve(v, boundaryConditions=BCs)
           >>> print numerix.allclose(v, answer, rtol=2e-5)
           True
        """
        pass

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
