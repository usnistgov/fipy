#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "ttri2Dinput.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:02:29 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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
`./examples/diffusion/steadyState/mesh1D/input.py`. The difference
being that it uses a triangular mesh loaded in using the GmshImporter.

The result is again tested in the same way:

    >>> Lx = 20
    >>> x = mesh.getCellCenters()[:,0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> import Numeric
    >>> Numeric.allclose(Numeric.array(var), analyticalArray, atol = 0.025)
    1

Note that this test case will only work if you run it by running the main FiPy
test suite. If you run it directly from the directory it is in it will not be able to find the mesh file.

"""

from fipy.meshes.grid2D import Grid2D
from fipy.equations.diffusionEquation import DiffusionEquation
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.iterators.iterator import Iterator
from fipy.variables.cellVariable import CellVariable
from fipy.viewers.pyxviewer import PyxViewer
from fipy.meshes.numMesh.gmshImport import GmshImporter2D
import sys
import Numeric

valueLeft = 0.
valueRight = 1.

mesh = GmshImporter2D("%s/%s" % (sys.__dict__['path'][0], "examples/diffusion/steadyState/mesh20x20/modifiedMesh.msh"))

var = CellVariable(name = "solution variable",
                   mesh = mesh,
                   value = valueLeft)

viewer = PyxViewer(var)

def leftSide(face):
    a = face.getCenter()[0]
    if(((a ** 2) < 0.000000000000001) and (face.getID() in mesh.getExteriorFaceIDs())):
        return 1
    else:
        return 0

def rightSide(face):
    a = face.getCenter()[0]
    if(( ((a - 20) ** 2) < 0.000000000000001) and (face.getID() in mesh.getExteriorFaceIDs())):
        return 1
    else:
        return 0

def bottomSide(face):
    a = face.getCenter()[1]
    if(((a ** 2) < 0.000000000000001) and (face.getID() in mesh.getExteriorFaceIDs())):
        return 1
    else:
        return 0

def topSide(face):
    a = face.getCenter()[0]
    if(( ((a - 20) ** 2) < 0.000000000000001) and (face.getID() in mesh.getExteriorFaceIDs())):
        return 1
    else:
        return 0

eq = DiffusionEquation(var,
                       transientCoeff = 0., 
                       diffusionCoeff = 1.,
                       solver = LinearPCGSolver(tolerance = 1.e-15, 
                                                steps = 1000
                                                ),
                       boundaryConditions = (FixedValue(mesh.getFacesWithFilter(leftSide), valueLeft),
                                             FixedValue(mesh.getFacesWithFilter(rightSide), valueRight),
                                             FixedFlux(mesh.getFacesWithFilter(topSide),0.),
                                             FixedFlux(mesh.getFacesWithFilter(bottomSide),0.)
                                             )
                       )

it = Iterator((eq,))

it.timestep()

if __name__ == '__main__':
    viewer.plot(resolution = 0.05)
    varArray = Numeric.array(var)
    x = mesh.getCellCenters()[:,0]
    analyticalArray = valueLeft + (valueRight - valueLeft) * x / 20
    errorArray = varArray - analyticalArray
    errorVar = CellVariable(name = "absolute error",
                   mesh = mesh,
                   value = abs(errorArray))
    errorViewer = PyxViewer(errorVar)
    errorViewer.plot(resolution = 0.05)
    NonOrthoVar = CellVariable(name = "non-orthogonality",
                               mesh = mesh,
                               value = mesh.getNonOrthogonality())
    NOViewer = PyxViewer(NonOrthoVar)    
    NOViewer.plot(resolution = 0.05)
    raw_input("finished")
