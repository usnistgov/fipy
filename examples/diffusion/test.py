"""Run all the test cases in examples/diffusion/
"""
from __future__ import unicode_literals

from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportDocTestSuite(testModuleNames = (
                                       'steadyState.test',
                                       'explicit.test',
                                       'nthOrder.test'
                                   ),
                                   docTestModuleNames = (
                                       'mesh1D',
                                       'mesh20x20',
                                       'circle',
                                       'circleQuad',
                                       'electrostatics',
                                       'variable',
                                       'anisotropy',
                                       'mesh20x20Coupled'
                                   ),
                                   base = __name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
