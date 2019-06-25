from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.viewers.matplotlibViewer.matplotlib2DViewer import AbstractMatplotlib2DViewer

__all__ = ["Matplotlib2DGridViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Matplotlib2DGridViewer(AbstractMatplotlib2DViewer):
    """
    Displays an image plot of a 2D `CellVariable` object using Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """

    __doc__ += AbstractMatplotlib2DViewer._test2D(viewer="Matplotlib2DGridViewer")

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
        float xmin, xmax, ymin, ymax, datamin, datamax : float, optional
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

        xmin, ymin = self.vars[0].mesh.extents['min']
        xmax, ymax = self.vars[0].mesh.extents['max']

        self.image = self.axes.imshow(self._data,
                                      extent=(xmin, xmax, ymin, ymax),
                                      vmin=0, vmax=1,
                                      cmap=self.cmap)

        self.axes.set_xlim(xmin=self._getLimit('xmin'),
                           xmax=self._getLimit('xmax'))

        self.axes.set_ylim(ymin=self._getLimit('ymin'),
                           ymax=self._getLimit('ymax'))

        if title is None:
            self.axes.set_title(self.vars[0].name)

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

    @property
    def _data(self):
        from fipy.tools.numerix import array, reshape
        return reshape(array(self.vars[0]), self.vars[0].mesh.shape[::-1])[::-1]

    def _plot(self):
        self.norm.vmin = self._getLimit(('datamin', 'zmin'))
        self.norm.vmax = self._getLimit(('datamax', 'zmax'))

        self.image.set_data(self.norm(self._data))

        if self.colorbar is not None:
            self.colorbar.plot()

def _test():
    from fipy.viewers.viewer import _test2D
    _test2D(Matplotlib2DGridViewer)

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
