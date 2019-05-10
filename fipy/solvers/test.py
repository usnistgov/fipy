from __future__ import unicode_literals
__all__ = []

from fipy.tests.lateImportTest import _LateImportTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportTestSuite(testModuleNames = (),
                                base = __name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
