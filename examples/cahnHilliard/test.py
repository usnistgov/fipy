from __future__ import unicode_literals
from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportDocTestSuite(docTestModuleNames = (
            'tanh1D',
            'mesh2D',
            'mesh3D',
            'sphere',
#             'mesh2DCoupled',  # FIXME: this test is [borked](https://github.com/usnistgov/fipy/issues/378)
        ), base = __name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
