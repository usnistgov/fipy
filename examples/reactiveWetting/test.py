"""Run all the test cases in `examples/reactiveWetting/`
"""
from __future__ import unicode_literals

from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportDocTestSuite(testModuleNames = (),
                                   docTestModuleNames = (
                                       'liquidVapor1D',
                                       'liquidVapor2D',
                                   ),
                                   base = __name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
