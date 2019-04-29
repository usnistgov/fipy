"""Run all the test cases in examples/meshing/
"""
from __future__ import unicode_literals

from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportDocTestSuite(testModuleNames = (),
                                   docTestModuleNames = (
                                       'sphere',
                                   ),
                                   base = __name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
