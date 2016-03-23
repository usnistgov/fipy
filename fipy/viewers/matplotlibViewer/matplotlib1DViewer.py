#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "matplotlib1DViewer.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.viewers.matplotlibViewer.matplotlibViewer import AbstractMatplotlibViewer

__all__ = ["Matplotlib1DViewer"]

class Matplotlib1DViewer(AbstractMatplotlibViewer):
    """
    Displays a y vs.  x plot of one or more 1D `CellVariable` objects using
    Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """

    __doc__ += AbstractMatplotlibViewer._test1D(viewer="Matplotlib1DViewer")

    def __init__(self, vars, title=None, xlog=False, ylog=False, limits={}, legend='upper left', axes=None, **kwlimits):
        """

        :Parameters:
          vars
            a `CellVariable` or tuple of `CellVariable` objects to plot
          title
            displayed at the top of the `Viewer` window
          xlog
            log scaling of x axis if `True`
          ylog
            log scaling of y axis if `True`
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          xmin, xmax, datamin, datamax
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
            (*ymin* and *ymax* are synonyms for *datamin* and *datamax*).
          legend
            place a legend at the specified position, if not `None`
          axes
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
        """
        kwlimits.update(limits)
        AbstractMatplotlibViewer.__init__(self, vars=vars, title=title, axes=axes, **kwlimits)

        import pylab

        if xlog and ylog:
            self.lines = [self.axes.loglog(*datum) for datum in self._data]
        elif xlog:
            self.lines = [self.axes.semilogx(*datum) for datum in self._data]
        elif ylog:
            self.lines = [self.axes.semilogy(*datum) for datum in self._data]
        else:
            self.lines = [self.axes.plot(*datum) for datum in self._data]

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

    def log():
        doc = "logarithmic data scaling"

        def fget(self):
            return self.axes.get_yscale() == 'log'

        def fset(self, value):
            ax = self.axes.get_yaxis()
            if value:
                ax = self.axes.set_yscale('log')
            else:
                ax = self.axes.set_yscale('linear')

        return locals()

    log = property(**log())

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

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
