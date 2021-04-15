r"""Solve the diffusion equation in a circular domain meshed with quadrangles.

This example demonstrates how to solve a simple diffusion problem on a
non-standard mesh with varying boundary conditions. The :term:`Gmsh` package
is used to create the mesh. Firstly, define some parameters for the
creation of the mesh,

>>> cellSize = 0.05
>>> radius = 1.

The `cellSize` is the preferred edge length of each mesh element and
the `radius` is the radius of the circular mesh domain. In the
following code section a file is created with the geometry that
describes the mesh. For details of how to write such geometry files
for :term:`Gmsh`, see the `gmsh manual`_.

.. _gmsh manual: http://www.geuz.org/gmsh/doc/texinfo/gmsh.html

The mesh created by :term:`Gmsh` is then imported into :term:`FiPy` using the
:class:`~fipy.meshes.gmshMesh.Gmsh2D` object.

>>> from fipy import CellVariable, Gmsh2D, TransientTerm, DiffusionTerm, Viewer
>>> from fipy.tools import numerix

>>> mesh = Gmsh2D('''
...               cellSize = %(cellSize)g;
...               radius = %(radius)g;
...               Point(1) = {0, 0, 0, cellSize};
...               Point(2) = {-radius, 0, 0, cellSize};
...               Point(3) = {0, radius, 0, cellSize};
...               Point(4) = {radius, 0, 0, cellSize};
...               Point(5) = {0, -radius, 0, cellSize};
...               Circle(6) = {2, 1, 3};
...               Circle(7) = {3, 1, 4};
...               Circle(8) = {4, 1, 5};
...               Circle(9) = {5, 1, 2};
...               Line Loop(10) = {6, 7, 8, 9};
...               Plane Surface(11) = {10};
...               Recombine Surface{11};
...               ''' % locals()) # doctest: +GMSH

Using this mesh, we can construct a solution variable

.. index::
   object: fipy.variables.cellVariable.CellVariable

>>> phi = CellVariable(name = "solution variable",
...                    mesh = mesh,
...                    value = 0.) # doctest: +GMSH

We can now create a :class:`Viewer <~fipy.viewers.viewer.AbstractViewer>` to see the mesh

>>> viewer = None
>>> from fipy import input
>>> if __name__ == '__main__':
...     try:
...         viewer = Viewer(vars=phi, datamin=-1, datamax=1.)
...         viewer.plotMesh()
...         input("Irregular circular mesh. Press <return> to proceed...")
...     except:
...         print("Unable to create a viewer for an irregular mesh (try Matplotlib2DViewer or MayaviViewer)")

.. image:: circleMesh.*
   :width: 90%
   :align: center

We set up a transient diffusion equation

.. index::
   object: fipy.terms.transientTerm.TransientTerm
   object: fipy.terms.implicitDiffusionTerm.DiffusionTerm

>>> D = 1.
>>> eq = TransientTerm() == DiffusionTerm(coeff=D)

The following line extracts the :math:`x` coordinate values on the exterior
faces. These are used as the boundary condition fixed values.

>>> X, Y = mesh.faceCenters # doctest: +GMSH
>>> phi.constrain(X, mesh.exteriorFaces) # doctest: +GMSH

We first step through the transient problem

>>> timeStepDuration = 10 * 0.9 * cellSize**2 / (2 * D)
>>> steps = 10
>>> from builtins import range
>>> for step in range(steps):
...     eq.solve(var=phi,
...              dt=timeStepDuration) # doctest: +GMSH
...     if viewer is not None:
...         viewer.plot()

.. image:: circleTransient.*
   :width: 90%
   :align: center

----

If we wanted to plot or analyze the results of this calculation with
another application, we could export tab-separated-values with

.. index::
   object: fipy.viewers.tsvViewer.TSVViewer

::

   TSVViewer(vars=(phi, phi.grad)).plot(filename="myTSV.tsv")

.. literalinclude:: myTSV.tsv

The values are listed at the :class:`~fipy.meshes.cell.Cell` centers.
Particularly for irregular meshes, no specific ordering should be relied upon.
Vector quantities are listed in multiple columns, one for each mesh dimension.

----

This problem again has an analytical solution that depends on the error
function, but it's a bit more complicated due to the varying boundary
conditions and the different horizontal diffusion length at different
vertical positions

>>> x, y = mesh.cellCenters # doctest: +GMSH
>>> t = timeStepDuration * steps

>>> phiAnalytical = CellVariable(name="analytical value",
...                              mesh=mesh) # doctest: +GMSH

.. index::
    module: scipy
    single: sqrt; arcsin; cos

>>> x0 = radius * numerix.cos(numerix.arcsin(y)) # doctest: +GMSH
>>> try:
...     from scipy.special import erf # doctest: +SCIPY
...     ## This function can sometimes throw nans on OS X
...     ## see http://projects.scipy.org/scipy/scipy/ticket/325
...     phiAnalytical.setValue(x0 * (erf((x0+x) / (2 * numerix.sqrt(D * t)))
...                                  - erf((x0-x) / (2 * numerix.sqrt(D * t))))) # doctest: +GMSH, +SCIPY
... except ImportError:
...     print("The SciPy library is not available to test the solution to \
... the transient diffusion equation")

>>> print(phi.allclose(phiAnalytical, atol = 7e-2)) # doctest: +GMSH, +SCIPY
1

>>> from fipy import input
>>> if __name__ == '__main__':
...     input("Transient diffusion. Press <return> to proceed...")

----

As in the earlier examples, we can also directly solve the steady-state
diffusion problem.

>>> DiffusionTerm(coeff=D).solve(var=phi) # doctest: +GMSH

The values at the elements should be equal to their `x` coordinate

>>> print(phi.allclose(x, atol = 0.035)) # doctest: +GMSH
1

Display the results if run as a script.

>>> from fipy import input
>>> if viewer is not None:
...     viewer.plot()
...     input("Steady-state diffusion. Press <return> to proceed...")

.. image:: circleSteadyState.*
   :width: 90%
   :align: center
"""
from __future__ import unicode_literals

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
