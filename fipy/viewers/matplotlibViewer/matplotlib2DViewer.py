#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "matplotlib2DViewer.py"
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

from fipy.viewers.matplotlibViewer.matplotlibViewer import AbstractMatplotlibViewer, _ColorBar

__all__ = ["Matplotlib2DViewer"]

class AbstractMatplotlib2DViewer(AbstractMatplotlibViewer):
    def figaspect(self, figaspect):
        if figaspect == 'auto':
            figaspect = self.vars[0].mesh.aspect2D
        return figaspect

class Matplotlib2DViewer(AbstractMatplotlib2DViewer):
    """
    Displays a contour plot of a 2D `CellVariable` object.

    The `Matplotlib2DViewer` plots a 2D `CellVariable` using Matplotlib_.

    .. _Matplotlib: http://matplotlib.sourceforge.net/
    """

    __doc__ += AbstractMatplotlib2DViewer._test2Dirregular(viewer="Matplotlib2DViewer")

    def __init__(self, vars, title=None, limits={}, cmap=None, colorbar='vertical', axes=None, figaspect='auto', **kwlimits):
        """Creates a `Matplotlib2DViewer`.


        :Parameters:
          vars
            a `CellVariable` object.
          title
            displayed at the top of the `Viewer` window
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          cmap
            the colormap. Defaults to `matplotlib.cm.jet`
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
          colorbar
            plot a colorbar in specified orientation if not `None`
          axes
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
          figaspect
            desired aspect ratio of figure. If arg is a number, use that aspect
            ratio. If arg is 'auto', the aspect ratio will be determined from
            the Variable's mesh.
        """
        kwlimits.update(limits)
        AbstractMatplotlib2DViewer.__init__(self, vars=vars, title=title, figaspect=figaspect,
                                            cmap=cmap, colorbar=colorbar, axes=axes,
                                            **kwlimits)

        self.mesh = self.vars[0].mesh

        vertexIDs = self.mesh._orderedCellVertexIDs

        vertexCoords = self.mesh.vertexCoords

        xCoords = numerix.take(vertexCoords[0], vertexIDs)
        yCoords = numerix.take(vertexCoords[1], vertexIDs)

        polys = []

        for x, y in zip(xCoords.swapaxes(0,1), yCoords.swapaxes(0,1)):
            if hasattr(x, 'mask'):
                x = x.compressed()
            if hasattr(y, 'mask'):
                y = y.compressed()
            polys.append(list(zip(x,y)))

        from matplotlib.collections import PolyCollection
        self.collection = PolyCollection(polys)
        self.collection.set_linewidth(0.5)
        try:
            self.axes.add_patch(self.collection)
        except:
            # PolyCollection not child of PatchCollection in matplotlib 0.98
            self.axes.add_collection(self.collection)

        xmin = self._getLimit('xmin', default=xCoords.min())
        xmax = self._getLimit('xmax', default=xCoords.max())
        ymin = self._getLimit('ymin', default=yCoords.min())
        ymax = self._getLimit('ymax', default=yCoords.max())

        self.axes.set_xlim(xmin=xmin, xmax=xmax)
        self.axes.set_ylim(ymin=ymin, ymax=ymax)

        self._plot()

    def _getSuitableVars(self, vars):
        from fipy.meshes.mesh2D import Mesh2D
        from fipy.variables.cellVariable import CellVariable
        vars = [var for var in AbstractMatplotlib2DViewer._getSuitableVars(self, vars) \
          if ((var.mesh.dim == 2 and isinstance(var, CellVariable))
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

        Z = self.vars[0].value

        self.norm.vmin = self._getLimit(('datamin', 'zmin'))
        self.norm.vmax = self._getLimit(('datamax', 'zmax'))

        rgba = self.cmap(self.norm(Z))

        self.collection.set_facecolors(rgba)
        self.collection.set_edgecolors(rgba)

        if self.colorbar is not None:
            self.colorbar.plot() #vmin=zmin, vmax=zmax)

##        pylab.xlim(xmin=self._getLimit('xmin'),
##                   xmax=self._getLimit('xmax'))

##        pylab.ylim(ymin=self._getLimit('ymin'),
##                   ymax=self._getLimit('ymax'))

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
