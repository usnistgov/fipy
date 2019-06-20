"""

This test file generates lots of different `SkewedGrid2D` meshes, each with a different non-orthogonality,
and runs a 1D diffusion problem on them all. It computes the RMS non-orthogonality and the RMS error
for each mesh and displays them in a graph, allowing the relationship of error to non-orthogonality to be investigated.
"""
from __future__ import division
from __future__ import unicode_literals

from builtins import range
if __name__ == '__main__':

    import sys
    import os

    from fipy import SkewedGrid2D, CellVariable, DiffusionTerm, Viewer
    from fipy.tools import numerix

    valueLeft = 0.
    valueRight = 1.

    meshList = []
    RMSNonOrthoList = []
    RMSErrorList = []

    for i in range(1, 501):
        meshList = meshList + [SkewedGrid2D(dx = 1.0, dy = 1.0, nx = 20, ny = 20, rand = (0.001 * i))]

    for mesh in meshList:
        var = CellVariable(name = "solution variable",
                           mesh = mesh,
                           value = valueLeft)

        var.constrain(valueLeft, mesh.facesLeft)
        var.constrain(valueRight, mesh.facesRight)

        DiffusionTerm().solve(var)

        varArray = numerix.array(var)
        x = mesh.cellCenters[0]
        analyticalArray = valueLeft + (valueRight - valueLeft) * x / 20
        errorArray = varArray - analyticalArray
        nonOrthoArray = mesh._nonOrthogonality
        RMSError = (numerix.add.reduce(errorArray * errorArray) / len(errorArray)) ** 0.5
        RMSNonOrtho = (numerix.add.reduce(nonOrthoArray * nonOrthoArray) / len(nonOrthoArray)) ** 0.5

        RMSNonOrthoList += [RMSNonOrtho]
        RMSErrorList += [RMSError]

    import pylab
    pylab.plot(RMSNonOrthoList, RMSErrorList, 'ro')
    pylab.show()
