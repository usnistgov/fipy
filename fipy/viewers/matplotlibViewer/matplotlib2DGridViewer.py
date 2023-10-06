from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.viewers.matplotlibViewer.abstractMatplotlib2DViewer import AbstractMatplotlib2DViewer

__all__ = ["Matplotlib2DGridViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Matplotlib2DGridViewer(AbstractMatplotlib2DViewer):
    """
    Displays an image plot of a 2D `CellVariable` object using Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """

    def __init__(self, vars, title=None, limits={}, cmap=None, colorbar='vertical', axes=None, figaspect='auto', **kwlimits):
        """Creates a `Matplotlib2DGridViewer`.

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable
            the `Variable` to display.
        title : str, optional
            displayed at the top of the `Viewer` window
        limits : dict, optional
            a (deprecated) alternative to limit keyword arguments
        cmap : ~matplotlib.colors.Colormap, optional
            the :class:`~matplotlib.colors.Colormap`.
            Defaults to `matplotlib.cm.jet`
        xmin, xmax, ymin, ymax, datamin, datamax : float, optional
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
        colorbar : bool, optional
            plot a color bar in specified orientation if not `None`
        axes : ~matplotlib.axes.Axes, optional
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
        figaspect : float, optional
            desired aspect ratio of figure. If arg is a number, use that aspect
            ratio. If arg is `auto`, the aspect ratio will be determined from
            the Variable's mesh.
        """
        kwlimits.update(limits)
        AbstractMatplotlib2DViewer.__init__(self, vars=vars, title=title,
                                            cmap=cmap, colorbar=colorbar, axes=axes, figaspect=figaspect,
                                            **kwlimits)

        self.axes.set_xlim(xmin=self._getLimit('xmin'),
                           xmax=self._getLimit('xmax'))

        self.axes.set_ylim(ymin=self._getLimit('ymin'),
                           ymax=self._getLimit('ymax'))

    def _getLimit(self, key, default=None):
        limit = AbstractMatplotlib2DViewer._getLimit(self, key, default=default)
        if limit is None:
            X, Y = self.vars[0].mesh.faceCenters
            if 'xmin' in key:
                limit = float(min(X))
            elif 'ymin' in key:
                limit = float(min(Y))
            if 'xmax' in key:
                limit = float(max(X))
            elif 'ymax' in key:
                limit = float(max(Y))
        return limit

    def _getSuitableVars(self, vars):
        from fipy.meshes.uniformGrid2D import UniformGrid2D
        from fipy.variables.cellVariable import CellVariable
        vars = [var for var in AbstractMatplotlib2DViewer._getSuitableVars(self, vars) \
          if (isinstance(var.mesh, UniformGrid2D) and isinstance(var, CellVariable)
              and var.rank == 0)]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError("Matplotlib2DGridViewer can only display a rank-0 CellVariable with a UniformGrid2D mesh")
        # this viewer can only display one variable
        return [vars[0]]

    def _make_mappable(self):
        xmin, ymin = self.vars[0].mesh.extents['min']
        xmax, ymax = self.vars[0].mesh.extents['max']

        image = self.axes.imshow(self._data,
                                 extent=(xmin, xmax, ymin, ymax),
                                 norm=self._norm,
                                 cmap=self.cmap)
        return image

    @property
    def _data(self):
        from fipy.tools.numerix import array, reshape
        return reshape(array(self.vars[0]), self.vars[0].mesh.shape[::-1])[::-1]

    def _plot(self):
        super(Matplotlib2DGridViewer, self)._plot()
        self._mappable.set_data(self._data)

    @classmethod
    def _doctest_body(cls):
        return cls._test2D()

def _test():
    from fipy.viewers.viewer import _test2D
    _test2D(Matplotlib2DGridViewer)

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
