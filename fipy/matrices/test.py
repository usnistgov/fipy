#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "test.py"
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

__all__ = []

from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram
from fipy.solvers import solver

if solver == 'trilinos':
    docTestModuleNames = ('trilinosMatrix', 'pysparseMatrix')
elif solver == 'no-pysparse':
    docTestModuleNames = ('trilinosMatrix',)
elif solver == 'scipy' or solver == 'pyamg':
    docTestModuleNames = ('scipyMatrix',)
elif solver == 'pysparse':
    docTestModuleNames = ('pysparseMatrix',)
else:
    raise ImportError('Unknown solver package %s' % solver)

def _suite():
    return _LateImportDocTestSuite(docTestModuleNames=docTestModuleNames, base=__name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
