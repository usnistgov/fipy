#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "vtkViewer.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Stiles  <daniel.stiles@nist.gov>
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


__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable

from fipy.viewers.vtkViewer.vtkViewer import VTKViewer

__all__ = ["VTKCellViewer"]

class VTKCellViewer(VTKViewer):
    """Renders `CellVariable` data in VTK format
    """
    def _makeDataSet(self, mesh):
        return mesh.VTKCellDataSet

    @property
    def _data(self):
        return self.dataset.cell_data

    @property
    def _variableClass(self):
        return CellVariable

    def _test(self):
        """
        >>> import os
        >>> from tempfile import mkstemp
        >>> f, fname = mkstemp(".vtk")
        >>> os.close(f)

        >>> try:
        ...     from tvtk.api import tvtk
        ... except ImportError, e:
        ...     from enthought.tvtk.api import tvtk
        ... # doctest: +TVTK

        >>> from fipy import *
        >>> from fipy.viewers.vtkViewer import VTKCellViewer

        >>> m = Grid1D(nx=10)
        >>> x, = m.cellCenters
        >>> v1 = CellVariable(mesh=m, value=x*x, name="x*x")
        >>> v2 = CellVariable(mesh=m, value=x)
        >>> v3 = v1.grad
        >>> v3.name = "v1.grad"
        >>> VTKCellViewer(vars=(v1, v2, v3)).plot(fname) # doctest: +TVTK
        >>> r = tvtk.DataSetReader() # doctest: +TVTK
        >>> r.file_name = fname # doctest: +TVTK
        >>> r.update() # doctest: +TVTK
        >>> c = r.output.cell_data # doctest: +TVTK
        >>> numerix.allclose(c.get_array("x*x").to_array(),
        ...                  v1.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(c.scalars.to_array(),
        ...                  v2.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(c.get_array("v1.grad").to_array().swapaxes(0,1)[0],
        ...                  v3.value) # doctest: +TVTK, +SERIAL
        True
        >>> r.get_vectors_name_in_file(0) == v3.name  # doctest: +TVTK, +PROCESSOR_0
        True

        >>> m = Grid2D(nx=1, ny=2)
        >>> x, y = m.cellCenters
        >>> v1 = CellVariable(mesh=m, value=x*y, name="x*y")
        >>> v2 = CellVariable(mesh=m, value=x*x) #, name="v2")
        >>> v3 = v1.grad
        >>> v3.name = "v1.grad"
        >>> VTKCellViewer(vars=(v1, v2, v3)).plot(fname) # doctest: +TVTK
        >>> r = tvtk.DataSetReader() # doctest: +TVTK
        >>> r.file_name = fname # doctest: +TVTK
        >>> r.update() # doctest: +TVTK
        >>> c = r.output.cell_data # doctest: +TVTK
        >>> numerix.allclose(c.get_array("x*y").to_array(),
        ...                  v1.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(c.scalars.to_array(),
        ...                  v2.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(c.get_array("v1.grad").to_array().swapaxes(0,1)[0:2],
        ...                  v3.value) # doctest: +TVTK, +SERIAL
        True
        >>> r.get_vectors_name_in_file(0) == v3.name  # doctest: +TVTK, +PROCESSOR_0
        True

        >>> m = (Grid2D(nx=5, ny=10, dx=0.1, dy=0.1)
        ...      + (Tri2D(nx=5, ny=5, dx=0.1, dy=0.1))
        ...      + ((0.5,), (0.2,)))
        >>> x, y = m.cellCenters
        >>> v1 = CellVariable(mesh=m, value=x*y, name="x*y")
        >>> v2 = CellVariable(mesh=m, value=x*x) #, name="v2")
        >>> v3 = v1.grad
        >>> v3.name = "v1.grad"
        >>> VTKCellViewer(vars=(v1, v2, v3)).plot(fname) # doctest: +TVTK
        >>> r = tvtk.DataSetReader() # doctest: +TVTK
        >>> r.file_name = fname # doctest: +TVTK
        >>> r.update() # doctest: +TVTK
        >>> c = r.output.cell_data # doctest: +TVTK
        >>> numerix.allclose(c.get_array("x*y").to_array(),
        ...                  v1.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(c.scalars.to_array(),
        ...                  v2.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(c.get_array("v1.grad").to_array().swapaxes(0,1)[0:2],
        ...                  v3.value) # doctest: +TVTK, +SERIAL
        True
        >>> r.get_vectors_name_in_file(0) == v3.name  # doctest: +TVTK, +PROCESSOR_0
        True

        >>> m = Grid3D(nx=2, ny=1, nz=1)
        >>> x, y, z = m.cellCenters
        >>> v1 = CellVariable(mesh=m, value=x*y*z, name="x*y*z")
        >>> v2 = CellVariable(mesh=m, value=x*y*y, name="x*y*y")
        >>> v3 = v1.grad
        >>> v3.name = "v1.grad"
        >>> VTKCellViewer(vars=(v1, v2, v3)).plot(filename=fname) # doctest: +TVTK
        >>> r = tvtk.DataSetReader() # doctest: +TVTK
        >>> r.file_name = fname # doctest: +TVTK
        >>> r.update() # doctest: +TVTK
        >>> c = r.output.cell_data # doctest: +TVTK
        >>> numerix.allclose(c.get_array("x*y*z").to_array(),
        ...                  v1.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(c.scalars.to_array(),
        ...                  v2.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(c.get_array("v1.grad").to_array().swapaxes(0,1),
        ...                  v3.value) # doctest: +TVTK, +SERIAL
        True
        >>> r.get_scalars_name_in_file(0) == v2.name  # doctest: +TVTK, +PROCESSOR_0
        True
        >>> r.get_vectors_name_in_file(0) == v3.name  # doctest: +TVTK, +PROCESSOR_0
        True

        >>> os.remove(fname)
        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
