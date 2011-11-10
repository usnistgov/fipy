#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mixedelement.py"
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
`./examples/diffusion/steadyState/mesh1D.py`. The difference in
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
    >>> steps = 10

A loop is required to execute the necessary time steps:

    >>> for step in range(steps):
    ...     eqn.solve(var, dt=timeStepDuration)
    
The result is again tested in the same way:

    >>> Lx = (2 * nx * dx)
    >>> x = bigMesh.cellCenters[0]
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    >>> ## print var.allclose(analyticalArray, rtol = 0.001, atol = 0.001)
    >>> print var.allclose(answer)
    1

"""
from fipy import *

dx = 1.
dy = 1.
nx = 10
ny = 2
valueLeft = 0.
valueRight = 1.
timeStepDuration = 0.005
steps = 10000

gridMesh = Grid2D(dx, dy, nx, ny)
## triMesh = Tri2D(dx, dy, nx, 1) + (dx*nx, 0)
triMesh = Tri2D(dx, dy, nx, 1) + ((dx*nx,), (0,))
bigMesh = gridMesh + triMesh

var = CellVariable(
    name = "concentration",
    mesh = bigMesh,
    value = valueLeft)

eqn = TransientTerm() == ExplicitDiffusionTerm()

exteriorFaces = bigMesh.exteriorFaces
xFace = bigMesh.faceCenters[0]

var.constrain(valueLeft, exteriorFaces & (xFace ** 2 < 0.000000000000001))
var.constrain((valueLeft + valueRight) * 0.5, exteriorFaces & ((xFace - (dx * nx)) ** 2 < 0.000000000000001))
var.constrain(valueRight, exteriorFaces & ((xFace - (2 * dx * nx)) ** 2 < 0.000000000000001))

answer = numerix.array([  0.00000000e+00,  8.78906250e-23,  1.54057617e-19,  1.19644866e-16,
        5.39556276e-14,  1.55308505e-11,  2.94461712e-09,  3.63798469e-07,
        2.74326174e-05,  1.01935828e-03,  9.76562500e-24,  1.92578125e-20,
        1.70937109e-17,  8.99433979e-15,  3.10726059e-12,  7.36603377e-10,
        1.21397338e-07,  1.37456643e-05,  1.02532568e-03,  4.57589878e-02,
        2.63278194e-07,  5.70863224e-12,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  1.51165440e-13,  1.23805218e-07,
        1.51873310e-03,  5.87457842e-01,  3.78270971e-06,  2.41898556e-10,
        2.62440000e-16,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  1.55453500e-09,  6.18653630e-05,  8.85109369e-02,
        7.24354518e-05,  1.32738123e-08,  8.11158300e-14,  0.00000000e+00,
        0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  2.27755930e-11,
        3.31776157e-06,  1.39387353e-02,  3.78270971e-06,  2.41898556e-10,
        2.62440000e-16,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        0.00000000e+00,  1.55453500e-09,  6.18653630e-05,  8.85109369e-02])

if __name__ == '__main__':
    viewer = Viewer(vars = var)
    for step in range(steps):
        var.updateOld()        
        eqn.solve(var, dt=timeStepDuration)
        if(not (step % 100)):
            print (step / 100)
    print var
    theMask = array([[10, 1, 20, 2]])
    viewer.plot()
##     viewer.plot(mask = theMask, graphwidth = 15, graphheight = 3)
    raw_input('finished')

