#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 10/7/04 {8:23:02 AM} { 5:14:21 PM}
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
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

r"""

This input example demonstrates how to create a non-standard mesh and
solve a simple diffusion example with varying boundary conditions. The
gmsh package is required to ruin this example. First set up some
parameters:

   >>> cellSize = 0.02
   >>> radius = 1.

The `cellSize` is the preffered edge length of each mesh element. The
meshing domain will be circular and thus a radius is defined. In the
following code section a file is created with the geometry that
describes the mesh. For details of how to write such geometry files
for gmsh, see the gmsh manual.

   >>> lines = [ 'cellSize = ' + str(cellSize) + ';\n',
   ...           'radius = ' + str(radius) + ';\n',
   ...           'Point(1) = {0, 0, 0, cellSize};\n',
   ...           'Point(2) = {-radius, 0, 0, cellSize};\n',
   ...           'Point(3) = {0, radius, 0, cellSize};\n',
   ...           'Point(4) = {radius, 0, 0, cellSize};\n',
   ...           'Point(5) = {0, -radius, 0, cellSize};\n',
   ...           'Circle(6) = {2, 1, 3};\n',
   ...           'Circle(7) = {3, 1, 4};\n',
   ...           'Circle(8) = {4, 1, 5};\n',
   ...           'Circle(9) = {5, 1, 2};\n',
   ...           'Line Loop(10) = {6, 7, 8, 9} ;\n',
   ...           'Plane Surface(11) = {10};\n']

   >>> import tempfile
   >>> (f, geomName) = tempfile.mkstemp('.geo')
   >>> file = open(geomName, 'w')
   >>> file.writelines(lines)
   >>> file.close()
   >>> import os
   >>> os.close(f)

The temporary file is used by gmsh to mesh the geometrically defined
region.

   >>> (f, meshName) = tempfile.mkstemp('.msh')
   >>> os.system('gmsh ' + geomName + ' -2 -v 0 -o ' + meshName)
   0
   >>> os.close(f)
   >>> os.remove(geomName)

The mesh created by gmsh is used to create a \FiPy{} mesh.
   
   >>> from fipy.meshes.numMesh.gmshImport import GmshImporter2D
   >>> mesh = GmshImporter2D(meshName)
   >>> os.remove(meshName)
    
Create a solution variable:

   >>> from fipy.variables.cellVariable import CellVariable
   >>> var = CellVariable(mesh = mesh, value = 0)

The following line extracts the x coordinate values on the exterior
faces. These are used as the boundary condition fixed values.

   >>> import Numeric
   >>> exteriorXcoords = Numeric.take(mesh.getFaceCenters()[:,0],
   ...                                   mesh.getExteriorFaceIDs())

The example is then solved as an implicit diffusion problem.

   >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
   >>> from fipy.boundaryConditions.fixedValue import FixedValue
   >>> ImplicitDiffusionTerm().solve(var = var,
   ...     boundaryConditions = (FixedValue(mesh.getExteriorFaces(), exteriorXcoords),))
                                                    
The values at the elements should be equal to the x coordinate

   >>> var.allclose(mesh.getCellCenters()[:,0], atol = 0.0079)
   1

Display the results if run as a script

   >>> if __name__ == '__main__':
   ...     from fipy.viewers.pyxviewer import PyxViewer
   ...     viewer = PyxViewer(var)
   ...     viewer.plot()

"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
    raw_input("finished")
