##from fipy.tools.profiler.profiler import Profiler
##from fipy.tools.profiler.profiler import calibrate_profiler

"""

This input file again solves a 1D diffusion problem as in
`./examples/diffusion/steadyState/mesh1D/input.py`. The difference
being that it uses a triangular mesh loaded in using the
:class:`~fipy.meshes.gmshMesh.Gmsh2D` object.

The result is again tested in the same way:

>>> from fipy import Gmsh2D, CellVariable, DiffusionTerm, Viewer
>>> from fipy.tools import numerix

>>> valueLeft = 0.
>>> valueRight = 1.

>>> Lx = 20
>>> mesh = Gmsh2D('''
...     cellSize = 0.5;
...     Point(2) = {0, 0, 0, cellSize};
...     Point(3) = {%(Lx)g, 0, 0, cellSize};
...     Point(4) = {%(Lx)g, %(Lx)g, 0, cellSize};
...     Point(5) = {0, %(Lx)g, 0, cellSize};
...
...     Line(6) = {2, 3};
...     Line(7) = {3, 4};
...     Line(8) = {4, 5};
...     Line(9) = {5, 2};
...
...     Line Loop(10) = {6, 7, 8, 9};
...
...     Plane Surface(11) = {10};
...     ''' % locals()) # doctest: +GMSH

>>> var = CellVariable(name = "solution variable",
...               mesh = mesh,
...               value = valueLeft) # doctest: +GMSH

>>> var.constrain(valueLeft, mesh.facesLeft) # doctest: +GMSH
>>> var.constrain(valueRight, mesh.facesRight) # doctest: +GMSH

>>> DiffusionTerm().solve(var) # doctest: +GMSH

>>> from fipy import input
>>> x = mesh.cellCenters[0] # doctest: +GMSH
>>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx # doctest: +GMSH
>>> print(var.allclose(analyticalArray, atol=0.025)) # doctest: +GMSH
True

>>> errorVar = abs(var - analyticalArray) # doctest: +GMSH
>>> errorVar.name = "absolute error" # doctest: +GMSH

>>> NonOrthoVar = CellVariable(name="non-orthogonality",
...                            mesh=mesh,
...                            value=mesh._nonOrthogonality) # doctest: +GMSH
>>> print(max(NonOrthoVar) < 0.51) # doctest: +GMSH
True

>>> if __name__ == '__main__':
...     viewer = Viewer(vars=var)
...     viewer.plot()
...
...     errorViewer = Viewer(vars=errorVar)
...     errorViewer.plot()
...
...     NOViewer = Viewer(vars=NonOrthoVar)
...     NOViewer.plot()
...
...     input("finished")
"""
from __future__ import unicode_literals


__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
