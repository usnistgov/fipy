from __future__ import unicode_literals
from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportDocTestSuite(docTestModuleNames = (
            'mesh1D.tri2Dinput',
            'mesh1D.inputPeriodic',
            'mesh20x20.tri2Dinput',
            'mesh50x50.input',
            'mesh50x50.tri2Dinput',
            'otherMeshes.grid3Dinput',
            'otherMeshes.prism',
            'mesh20x20.modifiedMeshInput',
            'mesh20x20.isotropy'
        ), base = __name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
