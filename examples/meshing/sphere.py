#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "sphere.py"
 #                                    created: 4/6/06 {11:26:11 AM}
 #                                last update: 10/5/07 {10:49:43 AM}
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://ctcms.nist.gov
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
 #  2006- 4- 6 JEG 1.0 original
 # ###################################################################
 ##

"""

An interesting problem is to solve an equation on a 2D geometry that
is embedded in 3D space, such as diffusion on the surface of a sphere
(with nothing either inside or outside the sphere). This example
demonstrates how to create the required mesh.

Test case.

   >>> max(numerix.sqrt(x**2 + y**2 + z**2)) < 5.3
   True
   >>> min(numerix.sqrt(x**2 + y**2 + z**2)) > 5.2
   True
   
"""

from fipy import *
import os

def dilate(x):
    return x * 1.1

mesh = GmshImporter2DIn3DSpace(os.path.splitext(__file__)[0] + '.gmsh').extrude(extrudeFunc=dilate)

x, y, z = mesh.getCellCenters()

var = CellVariable(mesh=mesh, value=x * y * z)

if __name__ == '__main__':

    viewer = MayaviViewer(vars = var, limits= {'datamin' : min(x * y * z), 'datamax' : max(x * y * z)})
    
    viewer.plot()

    raw_input('finished')
