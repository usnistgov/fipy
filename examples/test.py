


"""Run all the test cases in examples/
"""

from fipy.tests.lateImportTest import _LateImportTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportTestSuite(testModuleNames = (
        'diffusion.test',
        'chemotaxis.test',
        'phase.test',
        'convection.test',
        'elphf.test',
        'levelSet.test',
        'cahnHilliard.test',
        'flow.test',
        'meshing.test',
        'reactiveWetting.test',
        'riemann.test'
        ), base = __name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
