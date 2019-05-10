"""
Interactively test the viewers
"""
from __future__ import unicode_literals

from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram

__all__ = []

def _suite():
    return _LateImportDocTestSuite(testModuleNames = (
        'matplotlibViewer.test',
        'mayaviViewer.test',
        ), base = __name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
