#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "modifiedMeshInput.py"
 #
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
 # ###################################################################
 ##

"""

This input file again solves a 1D diffusion problem as in
`./examples/diffusion/steadyState/mesh1D/input.py`. The difference
being that it uses a triangular mesh loaded in using the Gmsh.

The result is again tested in the same way:

    >>> DiffusionTerm().solve(var, boundaryConditions = boundaryConditions)
    >>> Lx = 20
    >>> x = mesh.getCellCenters()[0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> print var.allclose(analyticalArray, atol = 0.026)
    1

Note that this test case will only work if you run it by running the
main FiPy test suite. If you run it directly from the directory it is
in it will not be able to find the mesh file.

"""

import sys

from fipy import *

valueLeft = 0.
valueRight = 1.

import os.path
mesh = Gmsh3D("""
    cellSize = 0.5;
    Len = 20;
    Hei = 1;
    Wid = 1;

    // first rectangle
    Point(1) = {0, 0, 0, cellSize};
    Point(2) = {0, 0, Wid, cellSize};
    Point(3) = {0, Hei, Wid, cellSize};
    Point(4) = {0, Hei, 0, cellSize};

    // second rectangle
    Point(5) = {Len, 0, 0, cellSize};
    Point(6) = {Len, 0, Wid, cellSize};
    Point(7) = {Len, Hei, Wid, cellSize};
    Point(8) = {Len, Hei, 0, cellSize};

    // lines joining first rect
    Line(9)  = {1, 2};
    Line(10) = {2, 3};
    Line(11) = {3, 4};
    Line(12) = {4, 1};

    // lines joining second rect
    Line(13) = {5, 6};
    Line(14) = {6, 7};
    Line(15) = {7, 8};
    Line(16) = {8, 5};

    // lines between the two "end caps"
    Line(17) = {1, 5};
    Line(18) = {2, 6};
    Line(19) = {3, 7};
    Line(20) = {4, 8};

    // lines around each side
    Line Loop(21) = {9, 10, 11, 12};
    Line Loop(22) = {13, 14, 15, 16};
    Line Loop(23) = {17, -16, -20, 12};
    Line Loop(24) = {13, -18, -9, 17};
    Line Loop(25) = {18, 14, -19, -10};
    Line Loop(26) = {-19, 11, 20, -15};

    Plane Surface(27) = {21};
    Plane Surface(28) = {22};
    Plane Surface(29) = {23};
    Plane Surface(30) = {24};
    Plane Surface(31) = {25};
    Plane Surface(32) = {26};

    Surface Loop(33) = {27, 28, 29, 30, 31, 32};

    Volume(34) = {33};
""")

var = CellVariable(name = "solution variable",
                   mesh = mesh,
                   value = valueLeft)

exteriorFaces = mesh.getExteriorFaces()
xFace = mesh.getFaceCenters()[0]
boundaryConditions = (FixedValue(exteriorFaces & (xFace ** 2 < 0.000000000000001), valueLeft),
                      FixedValue(exteriorFaces & ((xFace - 20) ** 2 < 0.000000000000001), valueRight))

if __name__ == '__main__':
    DiffusionTerm().solve(var, boundaryConditions = boundaryConditions)
    varArray = array(var)
    x = mesh.getCellCenters()[0]
    analyticalArray = valueLeft + (valueRight - valueLeft) * x / 20
    errorArray = varArray - analyticalArray
    errorVar = CellVariable(name = "absolute error",
                   mesh = mesh,
                   value = abs(errorArray))

