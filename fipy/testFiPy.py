"""Test suite for :term:`FiPy` modules
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'


__all__ = []

from fipy.tests.lateImportTest import _LateImportTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportTestSuite(testModuleNames = (
        'solvers.test',
        'terms.test',
        'tools.test',
        'matrices.test',
        'meshes.test',
        'variables.test',
        'viewers.test',
        'boundaryConditions.test',
    ), base = __name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
