#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "matplotlib2DContourViewer.py"
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

from fipy.tools import numerix

from fipy.viewers.matplotlibViewer.matplotlib2DViewer import AbstractMatplotlib2DViewer

__all__ = ["Matplotlib2DContourViewer"]

class Matplotlib2DContourViewer(AbstractMatplotlib2DViewer):
    """Displays a contour plot of a 2D `CellVariable` object.

    The `Matplotlib2DContourViewer` plots a 2D `CellVariable` using Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """

    __doc__ += Matplotlib2DContourViewer._test2D(viewer="Matplotlib2DContourViewer")


    def __init__(self, vars, title=None, limits={}, cmap=None, colorbar='vertical', axes=None, number=10, levels=None, figaspect='auto', **kwlimits):
        """Creates a `Matplotlib2DContourViewer`.

        :Parameters:
          vars
            a `CellVariable` object.
          title
            displayed at the top of the `Viewer` window
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
          cmap
            the colormap. Defaults to `matplotlib.cm.jet`
          colorbar
            plot a colorbar in specified orientation if not `None`
          axes
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
          number
            contour `number` automatically-chosen levels
          *levels* [level0, level1, ..., leveln]
            A list of floating point numbers indicating the level
            curves to draw; eg to draw just the zero contour pass
            ``levels=[0]``
          figaspect
            desired aspect ratio of figure. If arg is a number, use that aspect
            ratio. If arg is 'auto', the aspect ratio will be determined from
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
##         pylab.clf()

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
