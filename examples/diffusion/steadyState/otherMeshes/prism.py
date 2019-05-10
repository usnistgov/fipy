"""

This input file again solves a 1D diffusion problem as in
`./examples/diffusion/steadyState/mesh1D/input.py`. The difference
being that it uses a triangular mesh loaded in using the Gmsh.

The result is again tested in the same way:

    >>> from fipy import CellVariable, GmshGrid3D, DiffusionTerm

    >>> valueLeft = 0.
    >>> valueRight = 1.

    >>> mesh = GmshGrid3D(dx=1, dy=1, dz=1, nx=20, ny=1, nz=1) # doctest: +GMSH

    >>> var = CellVariable(name = "solution variable",
    ...                    mesh = mesh,
    ...                    value = valueLeft) # doctest: +GMSH

    >>> exteriorFaces = mesh.exteriorFaces # doctest: +GMSH
    >>> xFace = mesh.faceCenters[0] # doctest: +GMSH

    >>> var.constrain(valueLeft, exteriorFaces & (xFace ** 2 < 0.000000000000001)) # doctest: +GMSH
    >>> var.constrain(valueRight, exteriorFaces & ((xFace - 20) ** 2 < 0.000000000000001)) # doctest: +GMSH


    >>> DiffusionTerm().solve(var) # doctest: +GMSH
    >>> Lx = 20
    >>> x = mesh.cellCenters[0] # doctest: +GMSH
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx # doctest: +GMSH
    >>> print(var.allclose(analyticalArray, atol = 0.027)) # doctest: +GMSH
    1

"""
from __future__ import unicode_literals

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

