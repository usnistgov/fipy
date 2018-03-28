#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "orthoerror.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

"""

This test file generates lots of different SkewedGrid2D meshes, each with a different non-orthogonality,
and runs a 1D diffusion problem on them all. It ocmputes the RMS non-orthogonality and the RMS error
for each mesh and displays them in a graph, allowing the relationship of error to non-orthogonality to be investigated.
"""

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
