from __future__ import unicode_literals
from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportDocTestSuite(docTestModuleNames = (
            'exponential1D.mesh1D',
            'exponential1D.cylindricalMesh1D',
            'exponential1D.cylindricalMesh1DNonUniform',
            'exponential1D.tri2D',
            'exponential2D.mesh2D',
            'exponential2D.cylindricalMesh2D',
            'exponential2D.cylindricalMesh2DNonUniform',
            'exponential1DBack.mesh1D',
            'powerLaw1D.mesh1D',
            'exponential1DSource.mesh1D',
            'exponential2D.tri2D',
            'exponential1DSource.tri2D',
            'powerLaw1D.tri2D',
            'advection.vanLeerUpwind',
            'peclet',
            'robin',
            'source',
        ), base = __name__)

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
