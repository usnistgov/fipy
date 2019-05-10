"""
Test implementation of the viewers
"""
from __future__ import unicode_literals

__all__ = []

from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportDocTestSuite(testModuleNames = (
        'vtkViewer.test',),
                                   docTestModuleNames = (
        'tsvViewer',
        ), base = __name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
