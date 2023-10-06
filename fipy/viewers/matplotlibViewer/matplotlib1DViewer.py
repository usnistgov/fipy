from __future__ import unicode_literals
from builtins import zip
__docformat__ = 'restructuredtext'

from fipy.viewers.matplotlibViewer.abstractMatplotlibViewer import AbstractMatplotlibViewer

__all__ = ["Matplotlib1DViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Matplotlib1DViewer(AbstractMatplotlibViewer):
    """
    Displays a y vs.  x plot of one or more 1D `CellVariable` objects using
    Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """

    def __init__(self, vars, title=None, xlog=False, ylog=False, limits={}, legend='upper left', axes=None, **kwlimits):
        """

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable or list
            `CellVariable` objects to plot
        title : str, optional
            displayed at the top of the `Viewer` window
        xlog : bool
            log scaling of x axis if `True`
        ylog : bool
            log scaling of y axis if `True`
        limits : dict
            a (deprecated) alternative to limit keyword arguments
        xmin, xmax, datamin, datamax : float, optional
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
            (*ymin* and *ymax* are synonyms for *datamin* and *datamax*).
        legend : str
            place a legend at the specified position, if not `None`
        axes : ~matplotlib.axes.Axes
            if not `None`, `vars` will be plotted into this
            :ref:`Matplotlib` :class:`~matplotlib.axes.Axes` object
        """
        kwlimits.update(limits)
        AbstractMatplotlibViewer.__init__(self, vars=vars, title=title, axes=axes, **kwlimits)

        if xlog and ylog:
            self._lines = [self.axes.loglog(*datum) for datum in self._data]
        elif xlog:
            self._lines = [self.axes.semilogx(*datum) for datum in self._data]
        elif ylog:
            self._lines = [self.axes.semilogy(*datum) for datum in self._data]
        else:
            self._lines = [self.axes.plot(*datum) for datum in self._data]

        if legend is not None:
            self.axes.legend([var.name for var in self.vars], loc=legend)

        self.axes.set_xlim(xmin=self._getLimit('xmin'),
                           xmax=self._getLimit('xmax'))

        ymin = self._getLimit(('datamin', 'ymin'))
        ymax = self._getLimit(('datamax', 'ymax'))
        self.axes.set_ylim(ymin=ymin, ymax=ymax)

        if ymax is None or ymin is None:
            import warnings
            warnings.warn("Matplotlib1DViewer efficiency is improved by setting the 'datamax' and 'datamin' keys", UserWarning, stacklevel=2)

    @property
    def lines(self):
        """The collection of :ref:`Matplotlib` :class:`~matplotlib.lines.Line2D`
        objects representing the plotted data."""
        return self._lines

    @property
    def log(self):
        """Whether data has logarithmic scaling"""
        return self.axes.get_yscale() == 'log'

    @log.setter
    def log(self, value):
        ax = self.axes.get_yaxis()
        if value:
            ax = self.axes.set_yscale('log')
        else:
            ax = self.axes.set_yscale('linear')

    @property
    def _data(self):
        from fipy.tools.numerix import array
        return [[array(var.mesh.cellCenters[0]), array(var)] for var in self.vars]

    def _getSuitableVars(self, vars):
        vars = [var for var in AbstractMatplotlibViewer._getSuitableVars(self, vars) if var.mesh.dim == 1]

        if len(vars) > 1:
            vars = [var for var in vars if var.mesh is vars[0].mesh]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError("Can only plot 1D data")
        return vars

    def _plot(self):
        ymin, ymax = self._autoscale(vars=self.vars,
                                     datamin=self._getLimit(('datamin', 'ymin')),
                                     datamax=self._getLimit(('datamax', 'ymax')))

        self.axes.set_ylim(ymin=ymin, ymax=ymax)

        for line, datum in zip(self.lines, self._data):
            line[0].set_xdata(datum[0])
            line[0].set_ydata(datum[1])

    @classmethod
    def _doctest_body(cls):
        return cls._test1D()

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
