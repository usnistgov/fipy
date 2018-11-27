#!/usr/bin/env python


"""Run all the test cases in examples/
"""

from fipy.tests.lateImportTest import _LateImportTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportTestSuite(testModuleNames = (
            'distanceFunction.test',
            'advection.test',
            'surfactant.test',
            'electroChem.test',
        ), base = __name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
