from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = ["VTKViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

from fipy.viewers.viewer import AbstractViewer
from fipy.tests.doctestPlus import register_skipper

def _checkForTVTK():
    hasTVTK = True
    try:
        try:
            from tvtk.api import tvtk
        except ImportError as e:
            from enthought.tvtk.api import tvtk
    except Exception:
        hasTVTK = False
    return hasTVTK

register_skipper(flag="TVTK",
                 test=_checkForTVTK,
                 why="the `tvtk` package cannot be imported")

class VTKViewer(AbstractViewer):
    """Renders `_MeshVariable` data in VTK format
    """
    def __init__(self, vars, title=None, limits={}, **kwlimits):
        """Creates a `VTKViewer`

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable or ~fipy.variables.faceVariable.FaceVariable or list
            the `MeshVariable` objects to display.
        title : str, optional
            displayed at the top of the `Viewer` window
        limits : dict, optional
            a (deprecated) alternative to limit keyword arguments
        float xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax : float, optional
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
        """
        kwlimits.update(limits)
        AbstractViewer.__init__(self, vars=vars, title=title, **kwlimits)

        mesh = self.vars[0].mesh

        self.dataset = self._makeDataSet(mesh)

        data = self._data

        for var in self.vars:
            name, rank, value = self._nameRankValue(var)

            i = data.add_array(value)
            data.get_array(i).name = name

            if rank == 0:
                data.set_active_scalars(name)
            elif rank == 1:
                data.set_active_vectors(name)
            else:
                data.set_active_tensors(name)

    def _makeDataSet(self, mesh):
        pass

    @staticmethod
    def _nameRankValue(var):
        name = var.name or "%s #%d" % (var.__class__.__name__, id(var))
        rank = var.rank
        value = var.mesh._toVTK3D(var.value, rank=rank)

        return (name, rank, value)

    def plot(self, filename=None):
        data = self._data

        from fipy.tools import numerix

        for var in self.vars:
            name, rank, value = self._nameRankValue(var)

            if not (numerix.array(value.shape) == 0).any():
                data.get_array(name).to_array()[:] = value

        try:
            from tvtk.misc import write_data
        except ImportError as e:
            from enthought.tvtk.misc import write_data
        write_data(self.dataset, filename)

    def _getSuitableVars(self, vars):
        if type(vars) not in [type([]), type(())]:
            vars = [vars]
        cls = self._variableClass
        vars = [var for var in vars if isinstance(var, cls)]
        if len(vars) == 0:
            raise TypeError("%s can only display %s" % (self.__class__.__name__, cls.__name__))
        vars = [var for var in vars if var.mesh==vars[0].mesh]
        return vars


if __name__ == "__main__":
#     import fipy.tests.doctestPlus
#     fipy.tests.doctestPlus.execButNoTest()

    from fipy import *
    m = Grid3D(nx=2, ny=1, nz=1)
#     m = Grid3D(nx=3, ny=4, nz=5)
    x, y, z = m.cellCenters
    v1 = CellVariable(mesh=m, value=x*y*z, name="x*y*z")
    v2 = CellVariable(mesh=m, value=x*y*y, name="x*y*y")

    v3 = v1.grad
    v3.name = "v1.grad"
    v4 = v1.faceGrad
    v4.name = "v1.faceGrad"
    v5 = v1.harmonicFaceValue
    v5.name = "v1.harmonicFaceValue"
    v6 = v1.arithmeticFaceValue
    v6.name = "v1.arithmeticFaceValue"

#     vw = VTKViewer(vars=(v1, v2))
#     vw = VTKViewer(vars=(v1, v2, v3)) #, v4, v5, v6))
    vw = VTKViewer(vars=(v4, v5, v6))

    vw.plot(filename="face.vtk")

#     m = Grid2D(nx=1, ny=2)
#     x, y = m.cellCenters
#     v1 = CellVariable(mesh=m, value=x*y, name="v1")
#     v2 = CellVariable(mesh=m, value=x*x) #, name="v2")
#     vw = VTKViewer(vars=(v1, v2))

#     m = Grid1D(nx=10)
#     x,  = m.cellCenters
#     v1 = CellVariable(mesh=m, value=x*x, name="v1")
#     v2 = CellVariable(mesh=m, value=x) #, name="v2")
#     vw = VTKViewer(vars=(v1, v2))
