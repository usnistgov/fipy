#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "ttri2Dinput.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

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

>>> mesh = Gmsh2D('''
...     cellSize = 0.5;
...     Point(2) = {0, 0, 0, cellSize};
...     Point(3) = {20, 0, 0, cellSize};
...     Point(4) = {20, 20, 0, cellSize};
...     Point(5) = {0, 20, 0, cellSize};
...
...     Line(6) = {2, 3};
...     Line(7) = {3, 4};
...     Line(8) = {4, 5};
...     Line(9) = {5, 2};
...
...     Line Loop(10) = {6, 7, 8, 9};
...
...     Plane Surface(11) = {10};
...     ''') # doctest: +GMSH

>>> var = CellVariable(name = "solution variable",
...               mesh = mesh,
...               value = valueLeft) # doctest: +GMSH

>>> var.constrain(valueLeft, mesh.facesLeft) # doctest: +GMSH
>>> var.constrain(valueRight, mesh.facesRight) # doctest: +GMSH

>>> DiffusionTerm().solve(var) # doctest: +GMSH

>>> if __name__ == '__main__':
...     viewer = Viewer(vars = var)
...     viewer.plot()
...     varArray = numerix.array(var)
...     x = mesh.cellCenters[0]
...     analyticalArray = valueLeft + (valueRight - valueLeft) * x / 20
...     errorArray = varArray - analyticalArray
...     errorVar = CellVariable(name = "absolute error",
...                    mesh = mesh,
...                    value = abs(errorArray))
...     errorViewer = Viewer(vars = errorVar)
...     errorViewer.plot()
...
...     NonOrthoVar = CellVariable(name = "non-orthogonality",
...                           mesh = mesh,
...                           value = mesh._nonOrthogonality)
...     NOViewer = Viewer(vars = NonOrthoVar)
...
...
...     NOViewer.plot()
...     raw_input("finished")
... else:
...     Lx = 20
...     x = mesh.cellCenters[0] # doctest: +GMSH
...     analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx # doctest: +GMSH
...     print var.allclose(analyticalArray, atol = 0.025) # doctest: +GMSH
1

>>> max(mesh._nonOrthogonality) < 0.51 # doctest: +GMSH
True



Note that this test case will only work if you run it by running the
main FiPy test suite. If you run it directly from the directory it is
in it will not be able to find the mesh file.

"""


__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
