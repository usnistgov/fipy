from __future__ import division
from __future__ import unicode_literals
from builtins import zip
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.viewers.matplotlibViewer.abstractMatplotlibViewer import AbstractMatplotlibViewer
from fipy.viewers import MeshDimensionError

__all__ = ["MatplotlibSurfactantViewer"]

class MatplotlibSurfactantViewer(AbstractMatplotlibViewer):

    """
    The `MatplotlibSurfactantViewer` creates a viewer with the Matplotlib_ python
    plotting package that displays a `DistanceVariable`.

    .. _Matplotlib: http://matplotlib.sourceforge.net/

    """

    def __init__(self, distanceVar, surfactantVar=None, levelSetValue=0., title=None, smooth=0, zoomFactor=1., animate=False, limits={}, **kwlimits):
        """
        Create a `MatplotlibSurfactantViewer`.

            >>> from fipy import *
            >>> m = Grid2D(nx=100, ny=100)
            >>> x, y = m.cellCenters
            >>> v = CellVariable(mesh=m, value=x**2 + y**2 - 10**2)
            >>> s = CellVariable(mesh=m, value=sin(x / 10) * cos(y / 30))
            >>> viewer = MatplotlibSurfactantViewer(distanceVar=v, surfactantVar=s)
            >>> from builtins import range
            >>> for r in range(1, 200):
            ...     v.setValue(x**2 + y**2 - r**2)
            ...     viewer.plot()

            >>> from fipy import *
            >>> dx = 1.
            >>> dy = 1.
            >>> nx = 11
            >>> ny = 11
            >>> Lx = ny * dy
            >>> Ly = nx * dx
            >>> mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)
            >>> # from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
            >>> var = DistanceVariable(mesh = mesh, value = -1)

            >>> x, y = mesh.cellCenters

            >>> var.setValue(1, where=(x - Lx / 2.)**2 + (y - Ly / 2.)**2 < (Lx / 4.)**2)
            >>> var.calcDistanceFunction()
            >>> viewer = MatplotlibSurfactantViewer(var, smooth = 2)
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer

            >>> var = DistanceVariable(mesh = mesh, value = -1)

            >>> var.setValue(1, where=(y > 2. * Ly / 3.) | ((x > Lx / 2.) & (y > Ly / 3.)) | ((y < Ly / 6.) & (x > Lx / 2)))
            >>> var.calcDistanceFunction()
            >>> viewer = MatplotlibSurfactantViewer(var)
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer

            >>> viewer = MatplotlibSurfactantViewer(var, smooth = 2)
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer

        Parameters
        ----------
        distanceVar : ~fipy.variables.distanceVariable.DistanceVariable
        levelSetValue : float
            the value of the contour to be displayed
        title : str
            displayed at the top of the `Viewer` window
        animate : bool
            whether to show only the initial condition and the
            moving top boundary or to show all contours (Default)
        limits : dict, optional
            a (deprecated) alternative to limit keyword arguments
        xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax : float, optional
            displayed range of data. A 1D `Viewer` will only use `xmin` and
            `xmax`, a 2D viewer will also use `ymin` and `ymax`, and so on. All
            viewers will use `datamin` and `datamax`. Any limit set to a
            (default) value of `None` will autoscale.
        """

        kwlimits.update(limits)
        AbstractMatplotlibViewer.__init__(self, vars=[], title=title, **kwlimits)

        self.distanceVar = distanceVar
        if surfactantVar is None:
            self.surfactantVar = numerix.zeros(len(self.distanceVar), 'd')
        else:
            self.surfactantVar = surfactantVar
        self.smooth = smooth
        self.zoomFactor = zoomFactor

        self.animate = animate
        if animate:
            self._initialCondition = None

        if distanceVar.mesh.dim != 2:
            raise MeshDimensionError('The MatplotlibSurfactantViewer only works for 2D meshes.')

    def _plot(self):
        mesh = self.distanceVar.mesh
        shape = mesh.shape
        X, Y = mesh.cellCenters

        maxX = max(X)

        X = X.reshape(shape, order='F')
        Y = Y.reshape(shape, order='F')
        Z = self.distanceVar.value.reshape(shape, order='F')

        zmin, zmax = self._autoscale(vars=(self.surfactantVar,),
                                     datamin=self._getLimit(('datamin', 'zmin')),
                                     datamax=self._getLimit(('datamax', 'zmax')))

        import pylab
        import matplotlib

        CS = pylab.contour(X, Y, Z, (0.,))
        zc = CS.collections[0]

        verts = numerix.array(zc.get_verts())
        IDs = numerix.array([mesh._getNearestCellID(vert[..., numerix.newaxis]) for vert in verts])

        colors = pylab.cm.jet(( self.surfactantVar[IDs] - zmin) / (zmax - zmin))
        segments = list(zip(verts[:-1], verts[1:]))
        LC = matplotlib.collections.LineCollection(segments, colors=colors)

#         CS.ax.add_collection(LC)

        verts = numerix.array((-verts[..., 0], verts[..., 1])).swapaxes(0, 1)
        segments = list(zip(verts[:-1], verts[1:]))
        LC = matplotlib.collections.LineCollection(segments, colors=colors)

#         CS.ax.add_collection(LC)

        CS.ax.set_xlim((-maxX, maxX))

#         zc.set_color(pylab.cm.jet(IDs / mesh.numberOfCells))

#         pylab.xlim(xmin=self._getLimit('xmin'),
#                    xmax=self._getLimit('xmax'))
#
#         pylab.ylim(ymin=self._getLimit('ymin'),
#                    ymax=self._getLimit('ymax'))

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
