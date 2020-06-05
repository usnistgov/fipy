"""Test implementation of the mesh
"""
from __future__ import unicode_literals

__all__ = []

from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportDocTestSuite(docTestModuleNames = (
        'fipy.meshes.mesh',
        'fipy.meshes.mesh2D',
        'fipy.meshes.nonUniformGrid1D',
        'fipy.meshes.nonUniformGrid2D',
        'fipy.meshes.nonUniformGrid3D',
        'fipy.meshes.tri2D',
        'fipy.meshes.gmshMesh',
        'fipy.meshes.periodicGrid1D',
        'fipy.meshes.periodicGrid2D',
        'fipy.meshes.periodicGrid3D',
        'fipy.meshes.uniformGrid1D',
        'fipy.meshes.uniformGrid2D',
        'fipy.meshes.uniformGrid3D',
        'fipy.meshes.cylindricalUniformGrid1D',
        'fipy.meshes.cylindricalUniformGrid2D',
        'fipy.meshes.cylindricalNonUniformGrid1D',
        'fipy.meshes.cylindricalNonUniformGrid2D',
        'fipy.meshes.sphericalUniformGrid1D',
        'fipy.meshes.sphericalNonUniformGrid1D',
        'fipy.meshes.factoryMeshes',
        'fipy.meshes.abstractMesh',
        'fipy.meshes.representations.gridRepresentation'))

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
