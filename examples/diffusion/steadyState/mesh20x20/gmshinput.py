#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gmshinput.py"
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
`./examples/diffusion/steadyState/mesh1D/input.py`. In order to test the non-orthogonality error,
this uses a SkewedGrid2D, which is a Grid2D with each interior vertex moved in a random direction.
The RMS non-orthogonality and RMS error are measured for lots of different meshes and a graph is computed and
displayed. See the section of the manual under AdaptiveMesh for more information.

"""

from fipy.meshes.grid2D import Grid2D
from fipy.meshes.numMesh.skewedGrid2D import SkewedGrid2D
from fipy.meshes.numMesh.tri2D import Tri2D
from fipy.equations.diffusionEquation import DiffusionEquation
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.iterators.iterator import Iterator
from fipy.variables.cellVariable import CellVariable
from fipy.viewers.pyxviewer import PyxViewer
from fipy.meshes.numMesh.gmshImport import GmshImporter2D
import pyx
import sys
import os
import Numeric

valueLeft = 0.
valueRight = 1.

mesh = SkewedGrid2D(dx = 1.0, dy = 1.0, nx = 20, ny = 20, rand = 0.1)

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
    a = face.getCenter()[1]
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
                       boundaryConditions = (FixedValue(mesh.getFacesLeft(), valueLeft),
                                             FixedValue(mesh.getFacesRight(), valueRight),
                                             FixedFlux(mesh.getFacesTop(),0.),
                                             FixedFlux(mesh.getFacesBottom(),0.)
                                             )
                       )

it = Iterator((eq,))

it.timestep()

if __name__ == '__main__':
    ##viewer.plot(resolution = 0.1)
    varArray = Numeric.array(var)
    x = mesh.getCellCenters()[:,0]
    analyticalArray = valueLeft + (valueRight - valueLeft) * x / 20
    errorArray = varArray - analyticalArray
    errorVar = CellVariable(name = 'absolute error',
                            mesh = mesh,
                            value = abs(errorArray))
    errorViewer = PyxViewer(errorVar)
    ##errorViewer.plot(resolution = 0.1)
    NonOrthoVar = CellVariable(name = "non-orthogonality",
                               mesh = mesh,
                               value = mesh.getNonOrthogonality())
    NOViewer = PyxViewer(NonOrthoVar)    
    ##NOViewer.plot(resolution = 0.1)
    nonOrthoArray = mesh.getNonOrthogonality()
    displaylist = Numeric.concatenate(([nonOrthoArray], [errorArray]))
    displaylist = Numeric.transpose(displaylist)
    g = pyx.graph.graphxy(width = 8, x = pyx.graph.axis.linear(title = "Non-Orthogonality"), y = pyx.graph.axis.linear(title = "Error"))
    g.plot(pyx.graph.data.list(displaylist, addlinenumbers = 0, x=0, y=1))
    g.writeEPSfile("temptwo")
    os.system("gv temptwo.eps &")
    raw_input("finished")

