#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "linearGMRESSolver.py"
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

import sys

from fipy.solvers.scipy.scipySolver import ScipySolver

class LinearGMRESSolver(ScipySolver):
    """
    
    The `LinearGMRESSolver` solves a linear system of equations
    using the generalised minimal residual method (GMRES) with no
    GMRES solves systems with a general non-symmetric coefficient
    matrix.

    The `LinearGMRESSolver` is a wrapper class for the the
    Scipy_ `linalg.iterative.gmres()` method.

    .. warning::

        Currently the solvers that use Scipy_ are only useful for
        small systems due to the whole sparse matrix having to be
        turned into an array of size N * N.

    .. _Scipy: http://www.scipy.org

    
    """
    
    def _solve(self, L, x, b):
        """

        Tridiagonal test case,

           >>> N = 10
           >>> L = 1.
           >>> dx = L / N
           >>> from fipy.tools import numerix
           >>> a = numerix.zeros(N, 'd')
           >>> a[:] = 2 / dx
           >>> a[0] = 3 / dx
           >>> a[-1] = 3 / dx
           >>> solver = LinearGMRESSolver()
           >>> SparseMatrix = solver._getMatrixClass()
           >>> A = SparseMatrix(size = N)
           >>> A.addAtDiagonal(a)
           >>> ids = numerix.arange(N - 1)
           >>> A.addAt(-numerix.ones(N - 1, 'd') / dx, ids, ids + 1)
           >>> A.addAt(-numerix.ones(N - 1, 'd') / dx, ids + 1, ids)
           >>> b = numerix.zeros(N, 'd')
           >>> b[-1] = 2 / dx
           >>> x = numerix.zeros(N, 'd')
           >>> solver._solve(A, x, b)
           >>> numerix.allclose(x, numerix.arange(N) * dx + dx / 2.)
           1
           
        """
        from scipy.linalg.iterative import gmres
        x[:], info = gmres(L,b, x0 = x.copy(), tol = self.tolerance, maxiter = self.iterations)

        if (info != 0):
            print >> sys.stderr, 'gmres not converged'

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
