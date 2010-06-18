#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "circle.py"
 #
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

    >>> cellSize = 0.05
    >>> radius = 1.

The `cellSize` is the preferred edge length of each mesh element and
the `radius` is the radius of the circular mesh domain. In the
following code section a file is created with the geometry that
describes the mesh. For details of how to write such geometry files
for gmsh_, see the `gmsh manual`_.

.. _gmsh: http://www.geuz.org/gmsh/

.. _gmsh manual: http://www.geuz.org/gmsh/doc/texinfo/gmsh.html

The mesh created by gmsh_ is then imported into |FiPy| using the
`GmshImporter2D` object.
   
.. raw:: latex

   \IndexClass{GmshImporter2D}

..

    >>> from fipy import *
    >>> mesh = GmshImporter2D('''
    ...                       cellSize = %(cellSize)g;
    ...                       radius = %(radius)g;
    ...                       Point(1) = {0, 0, 0, cellSize};
    ...                       Point(2) = {-radius, 0, 0, cellSize};
    ...                       Point(3) = {0, radius, 0, cellSize};
    ...                       Point(4) = {radius, 0, 0, cellSize};
    ...                       Point(5) = {0, -radius, 0, cellSize};
    ...                       Circle(6) = {2, 1, 3};
    ...                       Circle(7) = {3, 1, 4};
    ...                       Circle(8) = {4, 1, 5};
    ...                       Circle(9) = {5, 1, 2};
    ...                       Line Loop(10) = {6, 7, 8, 9};
    ...                       Plane Surface(11) = {10};
    ...                       ''' % locals())
    
Using this mesh, we can construct a solution variable

.. raw:: latex

   \IndexClass{CellVariable}

..

    >>> phi = CellVariable(name = "solution variable",
    ...                    mesh = mesh,
    ...                    value = 0.)

We can now create a viewer to see the mesh

.. raw:: latex

   \IndexModule{viewers}

..

    >>> viewer = None
    >>> if __name__ == '__main__':
    ...     try:
    ...         viewer = Viewer(vars=phi, datamin=-1, datamax=1.)
    ...         viewer.plotMesh()
    ...         raw_input("Irregular circular mesh. Press <return> to proceed...")
    ...     except:
    ...         print "Unable to create a viewer for an irregular mesh (try Gist2DViewer or Matplotlib2DViewer)"

.. image:: examples/diffusion/circleMesh.pdf
   :scale: 50
   :align: center

We set up a transient diffusion equation

.. raw:: latex

   \IndexClass{TransientTerm}
   \IndexClass{ImplicitDiffusionTerm}

..

    >>> D = 1.
    >>> eq = TransientTerm() == ImplicitDiffusionTerm(coeff=D)

The following line extracts the `x` coordinate values on the exterior
faces. These are used as the boundary condition fixed values.

.. raw:: latex

   \IndexFunction{take}

..

    >>> X, Y = mesh.getFaceCenters()

.. raw:: latex

   \IndexClass{FixedValue}

..
    
    >>> BCs = (FixedValue(faces=mesh.getExteriorFaces(), value=X),)

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

    >>> x, y = mesh.getCellCenters()
    >>> t = timeStepDuration * steps

    >>> phiAnalytical = CellVariable(name="analytical value",
    ...                              mesh=mesh)

.. raw:: latex

   \IndexSoftware{SciPy}
   \IndexFunction{sqrt}
   \IndexFunction{arcsin}
   \IndexFunction{cos}

..

    >>> x0 = radius * cos(arcsin(y))
    >>> try:
    ...     from scipy.special import erf ## This function can sometimes throw nans on OS X
    ...                                   ## see http://projects.scipy.org/scipy/scipy/ticket/325
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

    >>> print phi.allclose(x, atol = 0.02)
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

