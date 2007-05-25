#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gmshinput.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 2/2/07 {8:51:09 AM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

"""

This input file again solves a 1D diffusion problem as in
`./examples/diffusion/steadyState/mesh1D/input.py`. In order to test
the non-orthogonality error, this uses a SkewedGrid2D, which is a
Grid2D with each interior vertex moved in a random direction.

"""

if __name__ == '__main__':
    import sys
    import os

    from fipy.tools import numerix

    from fipy.meshes.grid2D import Grid2D
    from fipy.meshes.skewedGrid2D import SkewedGrid2D
    from fipy.meshes.tri2D import Tri2D
    from fipy.boundaryConditions.fixedValue import FixedValue
    from fipy.boundaryConditions.fixedFlux import FixedFlux
    from fipy.variables.cellVariable import CellVariable
    import fipy.viewers
    from fipy.meshes.gmshImport import GmshImporter2D

    valueLeft = 0.
    valueRight = 1.

    mesh = SkewedGrid2D(dx = 1.0, dy = 1.0, nx = 20, ny = 20, rand = 0.1)

    var = CellVariable(name = "solution variable",
                       mesh = mesh,
                       value = valueLeft)

    viewer = fipy.viewers.make(vars = var)

    def leftSide(face):
        a = face.getCenter()[0]
        if(((a ** 2) < 0.000000000000001) and (face.getID() in mesh.getExteriorFaces())):
            return 1
        else:
            return 0

    def rightSide(face):
        a = face.getCenter()[0]
        if(( ((a - 20) ** 2) < 0.000000000000001) and (face.getID() in mesh.getExteriorFaces())):
            return 1
        else:
            return 0

    def bottomSide(face):
        a = face.getCenter()[1]
        if(((a ** 2) < 0.000000000000001) and (face.getID() in mesh.getExteriorFaces())):
            return 1
        else:
            return 0

    def topSide(face):
        a = face.getCenter()[1]
        if(( ((a - 20) ** 2) < 0.000000000000001) and (face.getID() in mesh.getExteriorFaces())):
            return 1
        else:
            return 0

    from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm

    ImplicitDiffusionTerm().solve(var, boundaryConditions = (FixedValue(mesh.getFacesLeft(), valueLeft),
                                                             FixedValue(mesh.getFacesRight(), valueRight)))

    varArray = numerix.array(var)
    x = mesh.getCellCenters()[0]
    analyticalArray = valueLeft + (valueRight - valueLeft) * x / 20
    errorArray = varArray - analyticalArray
    errorVar = CellVariable(name = 'absolute error',
                            mesh = mesh,
                            value = abs(errorArray))
    errorViewer = fipy.viewers.make(vars = errorVar)

    NonOrthoVar = CellVariable(name = "non-orthogonality",
                               mesh = mesh,
                               value = mesh._getNonOrthogonality())
    NOViewer = fipy.viewers.make(vars = NonOrthoVar)
    viewer.plot()
    NOViewer.plot()

    raw_input("finished")

