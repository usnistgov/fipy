__all__ = []

from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram
from fipy.solvers import solver_suite

docTestModuleNames = ['solver']
if solver_suite == 'scipy':
    docTestModuleNames.append('scipy.linearLUSolver')

def _suite():
    return _LateImportDocTestSuite(docTestModuleNames=docTestModuleNames,
                                   base=__name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
