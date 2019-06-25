from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.viewers.matplotlibViewer.matplotlib2DViewer import AbstractMatplotlib2DViewer

__all__ = ["Matplotlib2DContourViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Matplotlib2DContourViewer(AbstractMatplotlib2DViewer):
    """Displays a contour plot of a 2D `CellVariable` object.

    The `Matplotlib2DContourViewer` plots a 2D `CellVariable` using Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """

    __doc__ += AbstractMatplotlib2DViewer._test2D(viewer="Matplotlib2DContourViewer")


    def __init__(self, vars, title=None, limits={}, cmap=None, colorbar='vertical', axes=None, number=10, levels=None, figaspect='auto', **kwlimits):
        """Creates a `Matplotlib2DContourViewer`.

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable
            `Variable` to display
        title : str, optional
            displayed at the top of the `Viewer` window
        limits : dict, optional
            a (deprecated) alternative to limit keyword arguments
        float xmin, xmax, ymin, ymax, datamin, datamax : float, optional
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
        cmap : ~matplotlib.colors.Colormap, optional
            the Colormap.
            Defaults to `matplotlib.cm.jet`
        colorbar : bool, optional
            plot a color bar in specified orientation if not `None`
        axes : ~matplotlib.axes.Axes, optional
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
        number : int, optional
            contour `number` automatically-chosen levels
        levels : :obj:`list` of :obj:`float`, optional
            A list of numbers indicating the level
            curves to draw; e.g. to draw just the zero contour pass
            ``levels=[0]``
        figaspect : float
            desired aspect ratio of figure. If arg is a number, use that aspect
            ratio. If arg is `auto`, the aspect ratio will be determined from
            the Variable's mesh.

        """
        kwlimits.update(limits)
        AbstractMatplotlib2DViewer.__init__(self, vars=vars, title=title,
                                            cmap=cmap, colorbar=colorbar, axes=axes,
                                            figaspect=figaspect,
                                            **kwlimits)
        self.number = number
        self.levels = levels

        self._plot()

    def _getSuitableVars(self, vars):
        from fipy.meshes.mesh2D import Mesh2D
        from fipy.variables.cellVariable import CellVariable
        vars = [var for var in AbstractMatplotlib2DViewer._getSuitableVars(self, vars) \
          if ((isinstance(var.mesh, Mesh2D) and isinstance(var, CellVariable))
              and var.rank == 0)]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError("Matplotlib2DViewer can only display a rank-0, 2D CellVariable")
        # this viewer can only display one variable
        return [vars[0]]

    def _plot(self):
##         plt.clf()

##         ## Added garbage collection since matplotlib objects seem to hang
##         ## around and accumulate.
##         import gc
##         gc.collect()

        mesh = self.vars[0].mesh
        x, y = mesh.cellCenters
        z = self.vars[0].value

        xmin, ymin = mesh.extents['min']
        xmax, ymax = mesh.extents['max']

        from matplotlib.mlab import griddata

        xi = numerix.linspace(xmin, xmax, 1000)
        yi = numerix.linspace(ymin, ymax, 1000)
        # grid the data.
        zi = griddata(x, y, z, xi, yi, interp='linear')

        if hasattr(self, "_contourSet"):
            for collection in self._contourSet.collections:
                try:
                    ix = self.axes.collections.index(collection)
                except ValueError as e:
                    ix = None

                if ix is not None:
                    del self.axes.collections[ix]

        zmin, zmax = self._autoscale(vars=self.vars,
                                     datamin=self._getLimit(('datamin', 'zmin')),
                                     datamax=self._getLimit(('datamax', 'zmax')))

        self.norm.vmin = zmin
        self.norm.vmax = zmax

        if self.levels is not None:
            levels = self.levels
        else:
            levels = numerix.arange(self.number + 1) * (zmax - zmin) / self.number + zmin


        self._contourSet = self.axes.contour(xi, yi, zi, levels=levels, cmap=self.cmap)

        self.axes.set_xlim(xmin=self._getLimit('xmin'),
                           xmax=self._getLimit('xmax'))

        self.axes.set_ylim(ymin=self._getLimit('ymin'),
                           ymax=self._getLimit('ymax'))

        if self.colorbar is not None:
            self.colorbar.plot()


def _test():
    from fipy.viewers.viewer import _test2D
    _test2D(Matplotlib2DContourViewer)

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
