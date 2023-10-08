from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from future.builtins import super

from fipy.tools import numerix

from fipy.viewers.matplotlibViewer.abstractMatplotlib2DViewer import AbstractMatplotlib2DViewer

__all__ = ["Matplotlib2DGridContourViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Matplotlib2DGridContourViewer(AbstractMatplotlib2DViewer):
    """Displays a contour plot of a 2D `CellVariable` object.

    The `Matplotlib2DGridContourViewer` plots a 2D `CellVariable` using Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """

    def __init__(self, vars, title=None, limits={}, cmap=None, colorbar='vertical', axes=None, levels=None, figaspect='auto', **kwlimits):
        """Creates a `Matplotlib2DViewer`.

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable
            the `Variable` to display.
        title : str, optional
            displayed at the top of the `Viewer` window
        limits : dict
          a (deprecated) alternative to limit keyword arguments
        xmin, xmax, ymin, ymax, datamin, datamax : float, optional
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
        cmap : ~matplotlib.colors.Colormap, optional
            the :class:`~matplotlib.colors.Colormap`.
            Defaults to `matplotlib.cm.jet`
        colorbar : bool, optional
            plot a color bar if not `None`
        axes : ~matplotlib.axes.Axes, optional
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
        levels : int or array_like, optional
            Determines the number and positions of the contour lines /
            regions.  If an int `n`, tries to automatically choose no more
            than `n+1` "nice" contour levels over the range of `vars`.  If
            array_like, draw contour lines at the specified levels.  The
            values must be in increasing order.  E.g. to draw just the zero
            contour pass ``levels=[0]``.
        figaspect : float, optional
            desired aspect ratio of figure. If a number, use that aspect
            ratio. If `auto`, the aspect ratio will be determined from
            the *vars*'s mesh.
        """
        kwlimits.update(limits)
        AbstractMatplotlib2DViewer.__init__(self, vars=vars, title=title,
                                            cmap=cmap, colorbar=colorbar, axes=axes, figaspect=figaspect,
                                            **kwlimits)

        self.levels = levels

        self._plot()

    @property
    def levels(self):
        """The number of automatically-chosen contours or their values."""
        return self._levels

    @levels.setter
    def levels(self, value):
        self._levels = value

    def _getSuitableVars(self, vars):
        from fipy.meshes.nonUniformGrid2D import NonUniformGrid2D
        from fipy.meshes.uniformGrid2D import UniformGrid2D
        from fipy.variables.cellVariable import CellVariable
        vars = [var for var in AbstractMatplotlib2DViewer._getSuitableVars(self, vars) \
          if ((isinstance(var.mesh, NonUniformGrid2D)
               or isinstance(var.mesh, UniformGrid2D))
              and isinstance(var, CellVariable))]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError("The mesh must be a Grid2D instance")
        # this viewer can only display one variable
        return [vars[0]]

    def _plot(self):
        super(Matplotlib2DGridContourViewer, self)._plot()

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

        zmin = self._norm.vmin
        zmax = self._norm.vmax

        self.axes.contourf(X, Y, Z, levels=self.levels,
                           vmin=zmin, vmax=zmax, cmap=self.cmap)

        self.axes.set_xlim(xmin=self._getLimit('xmin'),
                           xmax=self._getLimit('xmax'))

        self.axes.set_ylim(ymin=self._getLimit('ymin'),
                           ymax=self._getLimit('ymax'))

    @classmethod
    def _doctest_body(cls):
        return cls._test2D()

    @classmethod
    def _doctest_extra(cls):
        return ("""
            >>> viewer.levels = 2
        """ + super()._doctest_extra())

def _test():
    from fipy.viewers.viewer import _test2D
    _test2D(Matplotlib2DGridContourViewer)

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
