#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "modifiedMeshInput.py"
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
    >>> print var.allclose(analyticalArray, atol = 0.027) # doctest: +GMSH
    1

"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
