#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "orthoerror.py"
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

"""

This test file generates lots of different SkewedGrid2D meshes, each with a different non-orthogonality,
and runs a 1D diffusion problem on them all. It ocmputes the RMS non-orthogonality and the RMS error
for each mesh and displays them in a graph, allowing the relationship of error to non-orthogonality to be investigated.
For more information, see the documentation for AdaptiveMesh.
"""

##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

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

meshList = []
displaylist = []
for i in range(1, 501):
    meshList = meshList + [SkewedGrid2D(dx = 1.0, dy = 1.0, nx = 20, ny = 20, rand = (0.001 * i))]


if __name__ == '__main__':
    
    for mesh in meshList:
        var = CellVariable(name = "solution variable",
                           mesh = mesh,
                           value = valueLeft)

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

        varArray = Numeric.array(var)
        x = mesh.getCellCenters()[:,0]
        analyticalArray = valueLeft + (valueRight - valueLeft) * x / 20
        errorArray = varArray - analyticalArray
        nonOrthoArray = mesh.getNonOrthogonality()
        RMSError = (Numeric.add.reduce(errorArray * errorArray) / len(errorArray)) ** 0.5
        RMSNonOrtho = (Numeric.add.reduce(nonOrthoArray * nonOrthoArray) / len(nonOrthoArray)) ** 0.5
        displaylist = displaylist + [[RMSNonOrtho, RMSError]]
    g = pyx.graph.graphxy(width = 8, xpos = 2, ypos = 2, x = pyx.graph.axis.linear(title = "RMS Non-Orthogonality"), y = pyx.graph.axis.linear(title = "RMS Error"))
    g.plot(pyx.graph.data.list(displaylist, addlinenumbers = 0, x=0, y=1), pyx.graph.style.symbol(size = "0.05 cm"))
    g.writeEPSfile("orthoerrorgraph")
    os.system("gv orthoerrorgraph.eps &")
    raw_input("finished")

