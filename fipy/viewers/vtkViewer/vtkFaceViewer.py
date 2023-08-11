from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.variables.faceVariable import FaceVariable

from fipy.viewers.vtkViewer.vtkViewer import VTKViewer

__all__ = ["VTKFaceViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class VTKFaceViewer(VTKViewer):
    """Renders :class:`~fipy.variables.meshVariable.MeshVariable` data in VTK format
    """
    def _makeDataSet(self, mesh):
        return mesh.VTKFaceDataSet

    @property
    def _data(self):
        return self.dataset.point_data

    @property
    def _variableClass(self):
        return FaceVariable

    def _test(self):
        """
        >>> import os
        >>> from tempfile import mkstemp
        >>> f, fname = mkstemp(".vtk")
        >>> os.close(f)

        >>> try:
        ...     from tvtk.api import tvtk
        ... except ImportError as e:
        ...     from enthought.tvtk.api import tvtk
        ... # doctest: +TVTK

        >>> from fipy import *
        >>> from fipy.viewers.vtkViewer import VTKFaceViewer

        >>> m = Grid1D(nx=10)
        >>> x, = m.cellCenters
        >>> v1 = CellVariable(mesh=m, value=x*x, name="x*x")
        >>> v2 = CellVariable(mesh=m, value=x)
        >>> v3 = v1.faceGrad
        >>> v3.name = "v1.faceGrad"
        >>> v4 = v1.harmonicFaceValue
        >>> v4.name = "v1.harmonicFaceValue"
        >>> v5 = v1.arithmeticFaceValue
        >>> v5.name = "v1.arithmeticFaceValue"
        >>> VTKFaceViewer(vars=(v3, v4, v5)).plot(fname) # doctest: +TVTK
        >>> r = tvtk.DataSetReader() # doctest: +TVTK
        >>> r.file_name = fname # doctest: +TVTK
        >>> r.update() # doctest: +TVTK
        >>> p = r.output.point_data # doctest: +TVTK
        >>> numerix.allclose(p.get_array("v1.faceGrad").to_array().swapaxes(0, 1)[0],
        ...                  v3.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(p.get_array("v1.harmonicFaceValue").to_array(),
        ...                  v4.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(p.get_array("v1.arithmeticFaceValue").to_array(),
        ...                  v5.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(p.scalars.to_array(),
        ...                  v5.value) # doctest: +TVTK, +SERIAL
        True
        >>> r.get_vectors_name_in_file(0) == v3.name  # doctest: +TVTK, +PROCESSOR_0
        True
        >>> r.get_scalars_name_in_file(0) == v5.name  # doctest: +TVTK, +PROCESSOR_0
        True

        >>> m = Grid2D(nx=1, ny=2)
        >>> x, y = m.cellCenters
        >>> v1 = CellVariable(mesh=m, value=x*y, name="x*y")
        >>> v2 = CellVariable(mesh=m, value=x*x) #, name="v2")
        >>> v3 = v1.faceGrad
        >>> v3.name = "v1.faceGrad"
        >>> v4 = v1.harmonicFaceValue
        >>> v4.name = "v1.harmonicFaceValue"
        >>> v5 = v1.arithmeticFaceValue
        >>> v5.name = "v1.arithmeticFaceValue"
        >>> VTKFaceViewer(vars=(v3, v4, v5)).plot(fname) # doctest: +TVTK
        >>> r = tvtk.DataSetReader() # doctest: +TVTK
        >>> r.file_name = fname # doctest: +TVTK
        >>> r.update() # doctest: +TVTK
        >>> p = r.output.point_data # doctest: +TVTK
        >>> numerix.allclose(p.get_array("v1.faceGrad").to_array().swapaxes(0, 1)[0:2],
        ...                  v3.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(p.get_array("v1.harmonicFaceValue").to_array(),
        ...                  v4.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(p.get_array("v1.arithmeticFaceValue").to_array(),
        ...                  v5.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(p.scalars.to_array(),
        ...                  v5.value) # doctest: +TVTK, +SERIAL
        True
        >>> r.get_vectors_name_in_file(0) == v3.name  # doctest: +TVTK, +PROCESSOR_0
        True
        >>> r.get_scalars_name_in_file(0) == v5.name  # doctest: +TVTK, +PROCESSOR_0
        True

        >>> m = (Grid2D(nx=5, ny=10, dx=0.1, dy=0.1)
        ...      + (Tri2D(nx=5, ny=5, dx=0.1, dy=0.1))
        ...      + ((0.5,), (0.2,)))
        >>> x, y = m.cellCenters
        >>> v1 = CellVariable(mesh=m, value=x*y, name="x*y")
        >>> v2 = CellVariable(mesh=m, value=x*x) #, name="v2")
        >>> v3 = v1.faceGrad
        >>> v3.name = "v1.faceGrad"
        >>> v4 = v1.harmonicFaceValue
        >>> v4.name = "v1.harmonicFaceValue"
        >>> v5 = v1.arithmeticFaceValue
        >>> v5.name = "v1.arithmeticFaceValue"
        >>> VTKFaceViewer(vars=(v3, v4, v5)).plot(fname) # doctest: +TVTK
        >>> r = tvtk.DataSetReader() # doctest: +TVTK
        >>> r.file_name = fname # doctest: +TVTK
        >>> r.update() # doctest: +TVTK
        >>> p = r.output.point_data # doctest: +TVTK
        >>> numerix.allclose(p.get_array("v1.faceGrad").to_array().swapaxes(0, 1)[0:2],
        ...                  v3.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(p.get_array("v1.harmonicFaceValue").to_array(),
        ...                  v4.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(p.get_array("v1.arithmeticFaceValue").to_array(),
        ...                  v5.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(p.scalars.to_array(),
        ...                  v5.value) # doctest: +TVTK, +SERIAL
        True
        >>> r.get_vectors_name_in_file(0) == v3.name  # doctest: +TVTK, +PROCESSOR_0
        True
        >>> r.get_scalars_name_in_file(0) == v5.name  # doctest: +TVTK, +PROCESSOR_0
        True

        >>> m = Grid3D(nx=2, ny=1, nz=1)
        >>> x, y, z = m.cellCenters
        >>> v1 = CellVariable(mesh=m, value=x*y*z, name="x*y*z")
        >>> v2 = CellVariable(mesh=m, value=x*y*y, name="x*y*y")
        >>> v3 = v1.faceGrad
        >>> v3.name = "v1.faceGrad"
        >>> v4 = v1.harmonicFaceValue
        >>> v4.name = "v1.harmonicFaceValue"
        >>> v5 = v1.arithmeticFaceValue
        >>> v5.name = "v1.arithmeticFaceValue"
        >>> VTKFaceViewer(vars=(v3, v4, v5)).plot(fname) # doctest: +TVTK
        >>> r = tvtk.DataSetReader() # doctest: +TVTK
        >>> r.file_name = fname # doctest: +TVTK
        >>> r.update() # doctest: +TVTK
        >>> p = r.output.point_data # doctest: +TVTK
        >>> numerix.allclose(p.get_array("v1.faceGrad").to_array().swapaxes(0, 1),
        ...                  v3.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(p.get_array("v1.harmonicFaceValue").to_array(),
        ...                  v4.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(p.get_array("v1.arithmeticFaceValue").to_array(),
        ...                  v5.value) # doctest: +TVTK, +SERIAL
        True
        >>> numerix.allclose(p.scalars.to_array(),
        ...                  v5.value) # doctest: +TVTK, +SERIAL
        True
        >>> r.get_vectors_name_in_file(0) == v3.name  # doctest: +TVTK, +PROCESSOR_0
        True
        >>> r.get_scalars_name_in_file(0) == v5.name  # doctest: +TVTK, +PROCESSOR_0
        True

        >>> os.remove(fname)
        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


