#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "circle.py"
 #                                    created: 4/6/06 {11:26:11 AM}
 #                                last update: 12/29/08 {10:31:53 PM}
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

r"""
This example demonstrates how to solve a simple diffusion problem on a
non-standard mesh with varying boundary conditions. The gmsh_ package
is used to create the mesh. Firstly, define some parameters for the
creation of the mesh,

.. raw:: latex

   \IndexSoftware{gmsh}
   
..

    >>> cellSize = 0.045
    >>> radius = 1.

The `cellSize` is the preferred edge length of each mesh element and
the `radius` is the radius of the circular mesh domain. In the
following code section a file is created with the geometry that
describes the mesh. For details of how to write such geometry files
for gmsh_, see the `gmsh manual`_.

.. _gmsh: http://www.geuz.org/gmsh/

.. _gmsh manual: http://www.geuz.org/gmsh/doc/texinfo/gmsh.html

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

The temporary file created above is used by gmsh_ to mesh the
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

The mesh created by gmsh_ is then imported into |FiPy| using the
`GmshImporter2D` object.
   
.. raw:: latex

   \IndexClass{GmshImporter2D}

..

    >>> from fipy.meshes.gmshImport import GmshImporter2D
    >>> mesh = GmshImporter2D(meshName)
    >>> os.remove(meshName)
    
Using this mesh, we can construct a solution variable

.. raw:: latex

   \IndexClass{CellVariable}

..

    >>> from fipy.variables.cellVariable import CellVariable
    >>> phi = CellVariable(name = "solution variable",
    ...                    mesh = mesh,
    ...                    value = 0)

We can now create a viewer to see the mesh (only the `Gist2DViewer` is
capable of displaying variables on this sort of irregular mesh)

.. raw:: latex

   \IndexSoftware{Pygist}
   \IndexSoftware{gist}
   \IndexClass{Gist2DViewer}
   \IndexModule{viewers}

..

    >>> viewer = None
    >>> if __name__ == '__main__':
    ...     try:
    ...         from fipy.viewers.gistViewer.gist2DViewer import Gist2DViewer
    ...         viewer = Gist2DViewer(vars=phi,
    ...                               limits={'datamin': -1, 'datamax': 1.})
    ...         viewer.plotMesh()
    ...         raw_input("Irregular circular mesh. Press <return> to proceed...")
    ...     except:
    ...         print "Unable to create a Gist2DViewer"

.. image:: examples/diffusion/circleMesh.pdf
   :scale: 50
   :align: center

We set up a transient diffusion equation

.. raw:: latex

   \IndexClass{TransientTerm}
   \IndexClass{ImplicitDiffusionTerm}

..

    >>> D = 1.
    >>> from fipy.terms.transientTerm import TransientTerm
    >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    >>> eq = TransientTerm() == ImplicitDiffusionTerm(coeff=D)

The following line extracts the `x` coordinate values on the exterior
faces. These are used as the boundary condition fixed values.

.. raw:: latex

   \IndexModule{numerix}
   \IndexFunction{take}

..

    >>> from fipy.tools import numerix
    >>> exteriorXcoords = numerix.take(mesh.getFaceCenters()[...,0],
    ...                                mesh.getExteriorFaces())

.. raw:: latex

   \IndexClass{FixedValue}

..
    
    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> BCs = (FixedValue(faces=mesh.getExteriorFaces(), value=exteriorXcoords),)

We first step through the transient problem

    >>> timeStepDuration = 10 * 0.9 * cellSize**2 / (2 * D)
    >>> steps = 10
    >>> for step in range(steps):
    ...     eq.solve(var=phi,
    ...              boundaryConditions=BCs,
    ...              dt=timeStepDuration)
    ...     if viewer is not None:
    ...         viewer.plot()

.. image:: examples/diffusion/circleTransient.pdf
   :scale: 50
   :align: center
   
-----

.. raw:: latex

   If we wanted to plot or analyze the results of this calculation with
   another application, we could export tab-separated-values with
   \IndexClass{TSVViewer}

::
    
   from fipy.viewers.tsvViewer import TSVViewer
   TSVViewer(vars=(phi, phi.getGrad())).plot(filename="myTSV.tsv")

.. raw:: latex

   \tabson[30] 
   {\tiny \verbatiminput{images/examples/diffusion/myTSV.tsv}} 
   \tabsoff
   
The values are listed at the `Cell` centers. Particularly for irregular
meshes, no specific ordering should be relied upon. Vector quantities are
listed in multiple columns, one for each mesh dimension.
            
-----

This problem again has an analytical solution that depends on the error
function, but it's a bit more complicated due to the varying boundary
conditions and the different horizontal diffusion length at different
vertical positions

    >>> x, y = mesh.getCellCenters()[...,0], mesh.getCellCenters()[...,1]
    >>> t = timeStepDuration * steps

    >>> phiAnalytical = CellVariable(name="analytical value",
    ...                              mesh=mesh)

.. raw:: latex

   \IndexModule{numerix}
   \IndexSoftware{SciPy}
   \IndexFunction{sqrt}
   \IndexFunction{arcsin}
   \IndexFunction{cos}

..

    >>> from fipy.tools.numerix import sqrt, arcsin, cos
    >>> x0 = radius * cos(arcsin(y))
    >>> try:
    ...     from scipy.special import erf
    ...     phiAnalytical.setValue(x0 * (erf((x0+x) / (2 * sqrt(D * t))) 
    ...                                  - erf((x0-x) / (2 * sqrt(D * t)))))
    ... except ImportError:
    ...     print "The SciPy library is not available to test the solution to \
    ... the transient diffusion equation"

    >>> print phi.allclose(phiAnalytical, atol = 7e-2)
    1

    >>> if __name__ == '__main__':
    ...     raw_input("Transient diffusion. Press <return> to proceed...")

-----

As in the earlier examples, we can also directly solve the steady-state
diffusion problem.

    >>> ImplicitDiffusionTerm(coeff=D).solve(var=phi,
    ...                                      boundaryConditions=BCs)
                                                    
The values at the elements should be equal to their `x` coordinate

    >>> print phi.allclose(mesh.getCellCenters()[...,0], atol = 0.02)
    1

Display the results if run as a script.

    >>> if viewer is not None:
    ...     viewer.plot()
    ...     raw_input("Steady-state diffusion. Press <return> to proceed...")

.. image:: examples/diffusion/circleSteadyState.pdf
   :scale: 50
   :align: center

.. |FiPy| raw:: latex

   \FiPy{}
"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

