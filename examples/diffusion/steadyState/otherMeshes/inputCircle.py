#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 6/7/05 {1:01:35 PM} { 5:14:21 PM}
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

This example demonstrates how to solve a simple diffusion problem on a
non-standard mesh with varying boundary conditions. The gmsh_ package
is used to create the mesh. Firstly, define some parameters for the
creation of the mesh,

   >>> cellSize = 0.05
   >>> radius = 1.

The `cellSize` is the preferred edge length of each mesh element and
the `radius` is the radius of the circular mesh domain. In the
following code section a file is created with the geometry that
describes the mesh. For details of how to write such geometry files
for gmsh_, see the gmsh_ manual_.

.. _gmsh: http://www.geuz.org/gmsh/

.. _manual: http://www.geuz.org/gmsh/doc/texinfo/gmsh.html

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

The temporary file created above is used by gmsh to mesh the
geometrically defined region.

   >>> import sys
   >>> if sys.platform == 'win32':
   ...     meshName = 'tmp.msh'
   ... else:
   ...     (f, meshName) = tempfile.mkstemp('.msh')
   >>> os.system('gmsh ' + geomName + ' -2 -v 0 -format msh -o ' + meshName)
   0
   >>> if sys.platform != 'win32':
   ...     os.close(f)
   >>> os.remove(geomName)

.. raw:: latex

   The mesh created by gmsh is then imported as a \FiPy{}

mesh using the `GmshImporter2D` object.
   
   >>> from fipy.meshes.gmshImport import GmshImporter2D
   >>> mesh = GmshImporter2D(meshName)
   >>> os.remove(meshName)
    
A solution variable is required.

   >>> from fipy.variables.cellVariable import CellVariable
   >>> var = CellVariable(mesh = mesh, value = 0)

The following line extracts the x coordinate values on the exterior
faces. These are used as the boundary condition fixed values.

   >>> import fipy.tools.array
   >>> exteriorXcoords = fipy.tools.array.take(mesh.getFaceCenters()[:,0],
   ...                                         mesh.getExteriorFaceIDs())

The example is then solved as an implicit diffusion problem.

   >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
   >>> from fipy.boundaryConditions.fixedValue import FixedValue
   >>> ImplicitDiffusionTerm().solve(var = var,
   ...     boundaryConditions = (
   ...         FixedValue(mesh.getExteriorFaces(), exteriorXcoords),))
                                                    
The values at the elements should be equal to the x coordinate

   >>> print var.allclose(mesh.getCellCenters()[:,0], atol = 0.02)
   1

Display the results if run as a script.

   >>> if __name__ == '__main__':
   ...     import fipy.viewers
   ...     fipy.viewers.make(vars = var).plotMesh()
   ...     fipy.viewers.make(vars = var,
   ...         limits = {'datamin': -1.0, 'datamax': 1.0}).plot()

.. image:: examples/diffusion/steadyState/otherMeshes/inputCircle.pdf
   :scale: 50
   :align: center

"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    raw_input("finished")
