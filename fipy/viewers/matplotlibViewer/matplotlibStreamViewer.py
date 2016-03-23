#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "matplotlibStreamViewer.py"
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
from fipy.variables.faceVariable import FaceVariable
from fipy.variables.cellVariable import CellVariable

from fipy.viewers.matplotlibViewer.matplotlib2DViewer import AbstractMatplotlib2DViewer

__all__ = ["MatplotlibStreamViewer"]

class MatplotlibStreamViewer(AbstractMatplotlib2DViewer):
    """Displays a stream plot of a 2D rank-1 `CellVariable` or
    `FaceVariable` object using Matplotlib_

    One issue is that this `Viewer` relies on `scipy.interpolate.griddata`,
    which interpolates on the convex hull of the data. The results is that
    streams are plotted across any concavities in the mesh.

    Another issue is that it does not seem possible to remove the streams
    without calling `cla()`, which means that different set of streams cannot be
    overlaid.

    .. _Matplotlib: http://matplotlib.sourceforge.net/

    """

    __doc__ += AbstractMatplotlib2DViewer._test2Dvector(viewer="MatplotlibStreamViewer")
    __doc__ += AbstractMatplotlib2DViewer._test2DvectorIrregular(viewer="MatplotlibStreamViewer")

    def __init__(self, vars, title=None, log=False, limits={}, axes=None, figaspect='auto',
                 density=1, linewidth=None, color=None, cmap=None, norm=None, arrowsize=1,
                 arrowstyle='-|>', minlength=0.1,
                 **kwlimits):
        """Creates a `MatplotlibStreamViewer`.

        :Parameters:
          vars
            a rank-1 `CellVariable` or `FaceVariable` object.
          title
            displayed at the top of the `Viewer` window
          log
            if `True`, arrow length goes at the base-10 logarithm of the magnitude
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
          axes
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
          figaspect
            desired aspect ratio of figure. If arg is a number, use that aspect
            ratio. If arg is 'auto', the aspect ratio will be determined from
            the Variable's mesh.
          *density* : float or 2-tuple
              Controls the closeness of streamlines. When `density = 1`, the domain
              is divided into a 25x25 grid---*density* linearly scales this grid.
              Each cell in the grid can have, at most, one traversing streamline.
              For different densities in each direction, use [density_x, density_y].
          linewidth : Numeric or rank-0 `MeshVariable`
            vary linewidth when given a `CellVariable` or `FaceVariable` of same
            type as vars.
          *color* : matplotlib color code, or rank-0 `MeshVariable`
              Streamline color. When given an array with the type as vars,
              *color* values are converted to colors using *cmap*.
          *cmap* : :class:`~matplotlib.colors.Colormap`
              Colormap used to plot streamlines and arrows. Only necessary when using
              an `MeshVariable` input for *color*.
          *norm* : :class:`~matplotlib.colors.Normalize`
              Normalize object used to scale luminance data to 0, 1. If None, stretch
              (min, max) to (0, 1). Only necessary when *color* is an `MeshVariable`.
          *arrowsize* : float
              Factor scale arrow size.
          *arrowstyle* : str
              Arrow style specification.
              See :class:`~matplotlib.patches.FancyArrowPatch`.
          *minlength* : float
              Minimum length of streamline in axes coordinates.

        """
        kwlimits.update(limits)
        AbstractMatplotlib2DViewer.__init__(self, vars=vars, title=title, axes=axes, figaspect=figaspect, **kwlimits)

        self.log = log
        self.kwargs = dict(density=density, cmap=cmap, norm=norm, arrowsize=arrowsize,
                           arrowstyle=arrowstyle, minlength=minlength)
        self.linewidth = linewidth
        self.color = color

        self._stream = None

        self._plot()

    def _getSuitableVars(self, vars):
        from fipy.meshes.mesh2D import Mesh2D
        from fipy.meshes.uniformGrid2D import UniformGrid2D

        vars = [var for var in AbstractMatplotlib2DViewer._getSuitableVars(self, vars) \
                if ((isinstance(var.mesh, Mesh2D)
                     or isinstance(var.mesh, UniformGrid2D))\
                    and (isinstance(var, FaceVariable) \
                         or isinstance(var, CellVariable)) and var.rank == 1)]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError("The mesh must be a Mesh2D instance")
        # this viewer can only display one variable
        return [vars[0]]

    def _plot(self):
        from scipy.interpolate import griddata

        var = self.vars[0]
        mesh = var.mesh

        xmin, ymin = mesh.extents['min']
        xmax, ymax = mesh.extents['max']

        N = 100
        X = numerix.linspace(xmin, xmax, N)
        Y = numerix.linspace(ymin, ymax, N)

        grid_x, grid_y = numerix.mgrid[xmin:xmax:N*1j, ymin:ymax:N*1j]

        if isinstance(var, FaceVariable):
            C = mesh.faceCenters
        elif isinstance(var, CellVariable):
            C = mesh.cellCenters

        U = griddata(C.value.T, var.value[0],
                     (grid_x, grid_y), method='cubic')
        V = griddata(C.value.T, var.value[1],
                     (grid_x, grid_y), method='cubic')

        lw = self.linewidth
        if isinstance(lw, (FaceVariable, CellVariable)):
            lw = griddata(C.value.T, lw.value,
                          (grid_x, grid_y), method='cubic')

        color = self.color
        if isinstance(color, (FaceVariable, CellVariable)):
            color = griddata(C.value.T, color.value,
                             (grid_x, grid_y), method='cubic', fill_value=color.min())

        U = U.T
        V = V.T

        ang = numerix.arctan2(V, U)
        mag = numerix.sqrt(U**2 + V**2)

        datamin, datamax = self._autoscale(vars=(mag,),
                                           datamin=self._getLimit('datamin'),
                                           datamax=self._getLimit('datamax'))

        mag = numerix.where(mag > datamax, numerix.nan, mag)
        mag = numerix.where(mag < datamin, numerix.nan, mag)

        if self.log:
            mag = numerix.log10(mag)

        U = mag * numerix.cos(ang)
        V = mag * numerix.sin(ang)

#         if self._stream is not None:
#             # the following doesn't work, nor does it help to `add_collection` first
#             # self._stream.arrows.remove()
#             self._stream.lines.remove()

        self.axes.cla()
        self._stream = self.axes.streamplot(X, Y, U, V, linewidth=lw, color=color, **self.kwargs)

        self.axes.set_xlim(xmin=self._getLimit('xmin'),
                           xmax=self._getLimit('xmax'))
        self.axes.set_ylim(ymin=self._getLimit('ymin'),
                           ymax=self._getLimit('ymax'))

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
