#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "pysparseSolver.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
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

import os

from fipy.matrices.scipyMatrix import _ScipyMeshMatrix
from fipy.solvers.pysparseMatrixSolver import _PysparseMatrixSolver

from fipy.tools.decorators import getsetDeprecated

class ScipySolver(_PysparseMatrixSolver):
    """
    The base `ScipySolver` class.
    
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, *args, **kwargs):
        if self.__class__ is ScipySolver:
            raise NotImplementedError, \
                  "can't instantiate abstract base class"
            
        super(ScipySolver, self).__init__(*args, **kwargs)
    
    @property
    def _matrixClass(self):
        return _ScipyMeshMatrix
                                   
    def _solve_(self, L, x, b):
        """
        Establishes a `pyamg.multilevel.multilevel_solver` object based on
        `self.solveFnc` and then solves, populating `relres` with
        a list of residuals.

        :Parameters:
            - `L`: A `fipy.matrices.pysparseMatrix._PysparseMeshMatrix`.
            - `x`: A `numpy.ndarray`.
            - `b`: A `numpy.ndarray`.
        """
        A = L.asformat('csr')
        M = None

        x, info = self.solveFnc(A, b, x, 
                                tol=self.tolerance, 
                                maxiter=self.iterations,
                                M=M)
        
        if os.environ.has_key('FIPY_VERBOSE_SOLVER'):
            from fipy.tools.debug import PRINT        
            PRINT('iterations: %d / %d' % (iter, self.iterations))
            
            if info < 0:
                PRINT('failure', self._warningList[info].__class__.__name__) 

        return x

