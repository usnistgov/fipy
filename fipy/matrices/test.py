from __future__ import unicode_literals
__all__ = []

from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram
from fipy.solvers import solver_suite

if solver_suite == 'trilinos':
    docTestModuleNames = ('trilinosMatrix', 'pysparseMatrix')
elif solver_suite == 'no-pysparse':
    docTestModuleNames = ('trilinosMatrix',)
elif solver_suite == 'scipy' or solver_suite == 'pyamg':
    docTestModuleNames = ('scipyMatrix',)
elif solver_suite == 'pysparse':
    docTestModuleNames = ('pysparseMatrix',)
elif solver_suite == 'pyamgx':
    docTestModuleNames = ('scipyMatrix',)
elif solver_suite == 'petsc':
    docTestModuleNames = ('petscMatrix',)
else:
    raise ImportError('Unknown solver package %s' % solver_suite)

def _suite():
    return _LateImportDocTestSuite(docTestModuleNames=docTestModuleNames, base=__name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
