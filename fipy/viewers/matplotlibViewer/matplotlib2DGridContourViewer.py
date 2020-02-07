from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.viewers.matplotlibViewer.matplotlib2DViewer import AbstractMatplotlib2DViewer

__all__ = ["Matplotlib2DGridContourViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Matplotlib2DGridContourViewer(AbstractMatplotlib2DViewer):
    """Displays a contour plot of a 2D `CellVariable` object.

    The `Matplotlib2DGridContourViewer` plots a 2D `CellVariable` using Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """

    __doc__ += AbstractMatplotlib2DViewer._test2D(viewer="Matplotlib2DGridContourViewer")


    def __init__(self, vars, title=None, limits={}, cmap=None, colorbar='vertical', axes=None, figaspect='auto', **kwlimits):
        """Creates a `Matplotlib2DViewer`.

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable
            the `Variable` to display.
        title : str, optional
            displayed at the top of the `Viewer` window
        limits : dict
          a (deprecated) alternative to limit keyword arguments
        float xmin, xmax, ymin, ymax, datamin, datamax : float, optional
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
        cmap : ~matplotlib.colors.Colormap, optional
            the :class:`~matplotlib.colors.Colormap`.
            Defaults to `matplotlib.cm.jet`
        colorbar : bool, optional
            plot a color bar if not `None`
        axes : ~matplotlib.axes.Axes, optional
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
        figaspect : float, optional
            desired aspect ratio of figure. If a number, use that aspect
            ratio. If `auto`, the aspect ratio will be determined from
            the *vars*'s mesh.
        """
        kwlimits.update(limits)
        AbstractMatplotlib2DViewer.__init__(self, vars=vars, title=title,
                                            cmap=cmap, colorbar=colorbar, axes=axes, figaspect=figaspect,
                                            **kwlimits)

        self._plot()

    def _getSuitableVars(self, vars):
        from fipy.meshes.grid2D import Grid2D
        from fipy.meshes.uniformGrid2D import UniformGrid2D
        from fipy.variables.cellVariable import CellVariable
        vars = [var for var in AbstractMatplotlib2DViewer._getSuitableVars(self, vars) \
          if ((isinstance(var.mesh, Grid2D)
               or isinstance(var.mesh, UniformGrid2D))
              and isinstance(var, CellVariable))]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError("The mesh must be a Grid2D instance")
        # this viewer can only display one variable
        return [vars[0]]

    def _plot(self):
##         plt.clf()

##         ## Added garbage collection since matplotlib objects seem to hang
##         ## around and accumulate.
##         import gc
##         gc.collect()

        mesh = self.vars[0].mesh
        shape = mesh.shape
        X, Y = mesh.cellCenters
        Z = self.vars[0].value
        X, Y, Z = [v.reshape(shape, order='F') for v in (X, Y, Z)]

        zmin, zmax = self._autoscale(vars=self.vars,
                                     datamin=self._getLimit(('datamin', 'zmin')),
                                     datamax=self._getLimit(('datamax', 'zmax')))

        self.norm.vmin = zmin
        self.norm.vmax = zmax

        numberOfContours = 10
        smallNumber = 1e-7
        diff = zmax - zmin

        if diff < smallNumber:
            V = numerix.arange(numberOfContours + 1) * smallNumber / numberOfContours + zmin
        else:
            V = numerix.arange(numberOfContours + 1) * diff / numberOfContours + zmin

        self.axes.contourf(X, Y, Z, V, cmap=self.cmap)

        self.axes.set_xlim(xmin=self._getLimit('xmin'),
                           xmax=self._getLimit('xmax'))

        self.axes.set_ylim(ymin=self._getLimit('ymin'),
                           ymax=self._getLimit('ymax'))

        if self.colorbar is not None:
            self.colorbar.plot()


def _test():
    from fipy.viewers.viewer import _test2D
    _test2D(Matplotlib2DGridContourViewer)

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
