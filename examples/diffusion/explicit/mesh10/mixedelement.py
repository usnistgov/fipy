#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:06:11 PM} 
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

"""

This input file again solves a 1D diffusion problem as in
`./examples/diffusion/steadyState/mesh1D/input.py`. The difference in
this example is that the solution method is explicit. The equation
used is the `ExplicitDiffusionEquation`. In this case many steps have
to be taken to reach equilibrum. The `timeStepDuration` parameter
specifies the size of each time step and `steps` is the number of
time steps.

    >>> dx = 1.
    >>> dy = 1.
    >>> nx = 10
    >>> ny = 2
    >>> valueLeft = 0.
    >>> valueRight = 1.
    >>> timeStepDuration = 0.005
    >>> steps = 25000

A loop is required to execute the necessary time steps:

    >>> for step in range(steps):
    ...     it.timestep()
    
The result is again tested in the same way:

    >>> Lx = (2 * nx * dx)
    >>> x = bigMesh.getCellCenters()[:,0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> import Numeric
    >>> Numeric.allclose(Numeric.array(var), analyticalArray, rtol = 0.001, atol = 0.001)
    1

"""



from fipy.meshes.grid2D import Grid2D
from fipy.meshes.numMesh.tri2D import Tri2D
from fipy.equations.explicitDiffusionEquation import ExplicitDiffusionEquation
from fipy.solvers.linearPCGSolver import LinearPCGSolver
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.iterators.iterator import Iterator
from fipy.variables.cellVariable import CellVariable
from fipy.viewers.pyxviewer import PyxViewer
import Numeric

dx = 1.
dy = 1.
nx = 10
ny = 2
valueLeft = 0.
valueRight = 1.
timeStepDuration = 0.005
steps = 25000

gridMesh = Grid2D(dx, dy, nx, ny)
triMesh = Tri2D(dx, dy, nx, 1) + (dx*nx, 0)
bigMesh = gridMesh + triMesh

## filter functions

def leftSide(face):
    a = face.getCenter()[0]
    if(((a ** 2) < 0.000000000000001) and (face.getID() in bigMesh.getExteriorFaceIDs())):
        return 1
    else:
        return 0

def inMiddle(face):
    a = face.getCenter()[0]
    if(( ((a - (dx * nx)) ** 2) < 0.000000000000001) and (face.getID() in bigMesh.getExteriorFaceIDs())):
        return 1
    else:
        return 0

def rightSide(face):
    a = face.getCenter()[0]
    if(( ((a - (2 * dx * nx)) ** 2) < 0.000000000000001) and (face.getID() in bigMesh.getExteriorFaceIDs())):
        return 1
    else:
        return 0

def allOthers(face):
    
    if((leftSide(face) or inMiddle(face) or rightSide(face)) or not (face.getID() in bigMesh.getExteriorFaceIDs())):
        return 0
    else:
        return 1   

var = CellVariable(
    name = "concentration",
    mesh = bigMesh,
    value = valueLeft)


eq = ExplicitDiffusionEquation(
    var,
    transientCoeff = 1. / timeStepDuration, 
    diffusionCoeff = 1.,
    solver = LinearPCGSolver(
    tolerance = 1.e-15, 
    steps = 1000
    ),
    boundaryConditions=(
    FixedValue(bigMesh.getFacesWithFilter(leftSide), valueLeft),
    FixedValue(bigMesh.getFacesWithFilter(inMiddle), (valueLeft + valueRight) * 0.5),
    FixedValue(bigMesh.getFacesWithFilter(rightSide), valueRight),
    FixedFlux(bigMesh.getFacesWithFilter(allOthers), 0.)
    )
    )

it = Iterator((eq,))

if __name__ == '__main__':
    viewer = PyxViewer(var)
    for step in range(steps):
        it.timestep()
        if(not (step % 100)):
            print (step / 100)
    theMask = Numeric.array([[10, 1, 20, 2]])
    viewer.plot(mask = theMask, graphwidth = 15, graphheight = 3)
    raw_input('finished')

