"""Test numeric implementation of the mesh
"""
from __future__ import unicode_literals

__all__ = []

from fipy.tests.doctestPlus import _LateImportDocTestSuite
import fipy.tests.testProgram

def _suite():
    return _LateImportDocTestSuite(
        docTestModuleNames = (
            'fipy.variables.variable',
            'fipy.variables.meshVariable',
            'fipy.variables.cellVariable',
            'fipy.variables.faceVariable',
            'fipy.variables.operatorVariable',
            'fipy.variables.betaNoiseVariable',
            'fipy.variables.exponentialNoiseVariable',
            'fipy.variables.gammaNoiseVariable',
            'fipy.variables.gaussianNoiseVariable',
            'fipy.variables.uniformNoiseVariable',
            'fipy.variables.modularVariable',
            'fipy.variables.binaryOperatorVariable',
            'fipy.variables.unaryOperatorVariable',
            'fipy.variables.coupledCellVariable',
            'fipy.variables.cellToFaceVariable',
            'fipy.variables.faceGradVariable',
            'fipy.variables.gaussCellGradVariable',
            'fipy.variables.faceGradContributionsVariable',
            'fipy.variables.surfactantConvectionVariable',
            'fipy.variables.surfactantVariable',
            'fipy.variables.levelSetDiffusionVariable',
            'fipy.variables.distanceVariable'
        ))

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')
