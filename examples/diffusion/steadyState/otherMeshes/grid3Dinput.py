#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "ttri2Dinput.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 7/5/07 {6:33:59 PM} 
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
Test case for the Grid3D.
 
   >>> ImplicitDiffusionTerm().solve(var, boundaryConditions = boundaryConditions)
   >>> ImplicitDiffusionTerm().solve(var2, boundaryConditions = boundaryConditions2)
   >>> a = numerix.array(var)
   >>> b = numerix.array(var2)
   >>> c = numerix.ravel(numerix.array((b, b, b)))
   >>> print numerix.allclose(a, c)
   1
   
"""

from fipy import *

nx = 10
ny = 5
nz = 3

dx = 1.
dy = 1.
dz = 1.

valueBottomTop = 0.
valueLeftRight = 1.

mesh = Grid3D(dx = dx, dy = dy, dz = dz, nx = nx, ny = ny, nz = nz)

var = CellVariable(name = "solution variable",
                   mesh = mesh,
                   value = valueBottomTop)

boundaryConditions = (FixedValue(mesh.getFacesLeft(),valueLeftRight),
                      FixedValue(mesh.getFacesRight(),valueLeftRight),
                      FixedValue(mesh.getFacesTop(),valueBottomTop),
                      FixedValue(mesh.getFacesBottom(),valueBottomTop))
                      
#do the 2D problem for comparison

nx = 10
ny = 5

dx = 1.
dy = 1.

mesh2 = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)

var2 = CellVariable(name = "solution variable 2D",
                    mesh = mesh2,
                    value = valueBottomTop)

boundaryConditions2 = (FixedValue(mesh2.getFacesLeft(),valueLeftRight),
                       FixedValue(mesh2.getFacesRight(),valueLeftRight),
                       FixedValue(mesh2.getFacesTop(),valueBottomTop),
                       FixedValue(mesh2.getFacesBottom(),valueBottomTop))

eqn = ImplicitDiffusionTerm()  

if __name__ == '__main__':
    eqn.solve(var2, boundaryConditions = boundaryConditions2)
    viewer = viewers.make(var2)
    viewer.plot()
    raw_input("finished")
