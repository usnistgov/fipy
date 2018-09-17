#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "inputGrid2D.py"
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

r"""

To run this example from the base FiPy directory, type::

    $ python examples/meshing/inputGrid2D.py --numberOfElements=X

This example demonstrates how to build a 1D mesh and obtain basic mesh
information. The command line argument, X, controls the number of
elements on the mesh. Firstly parse the command line argument for
`numberOfElements`, with the default set at 100.

   >>> from fipy.tools.parser import parse
   >>> numberOfElements = parse('--numberOfElements', action = 'store', type = 'int', default = 100)

A `Grid2D` object is invoked in the following way,

   >>> from fipy import Grid2D, CellVariable, Viewer
   >>> from fipy.tools import numerix

   >>> nx = int(numerix.sqrt(numberOfElements))
   >>> ny = nx
   >>> dx = 1.
   >>> dy = 1.
   >>> mesh = Grid2D(nx = nx, ny = nx, dx = dx, dy = dy)

Once the mesh has been built information about the mesh can be
obtained.  For example the mesh volumes can be obtained with the
`getCellVolumes()` method.

   >>> vols = mesh.cellVolumes
   >>> numerix.allclose(dx * dy * numerix.ones(nx * ny), vols)
   1

Obtain the number of cells in the mesh

   >>> N = mesh.numberOfCells
   >>> numerix.allclose(N, numberOfElements)
   1

Obtain all the left exterior faces, this is equal to `ny`.

   >>> faces = mesh.facesLeft
   >>> len(faces) == ny
   1

One can view the mesh with the following code,

   >>> if __name__ == '__main__':
   ...     viewer = Viewer(CellVariable(value = 0, mesh = mesh))
   ...     viewer.plot()

"""
__docformat__ = 'restructuredtext'

def _run():
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript(__name__))

if __name__ == '__main__':
    _run()
    raw_input("finished")
