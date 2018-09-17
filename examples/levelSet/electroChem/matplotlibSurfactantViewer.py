#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "mayaviSurfactantViewer.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 #  Description:
 #
 #  History
 #
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.viewers.matplotlibViewer.matplotlibViewer import AbstractMatplotlibViewer
from fipy.viewers import MeshDimensionError

__all__ = ["MatplotlibSurfactantViewer"]

class MatplotlibSurfactantViewer(AbstractMatplotlibViewer):

    """
    The `MatplotlibSurfactantViewer` creates a viewer with the Matplotlib_ python
    plotting package that displays a `DistanceVariable`.

    .. _Matplotlib: http://matplotlib.sourceforge.net/

    """

    def __init__(self, distanceVar, surfactantVar=None, levelSetValue=0., title=None, smooth=0, zoomFactor=1., animate=False, limits={}, **kwlimits):
        """
        Create a `MatplotlibSurfactantViewer`.

            >>> from fipy import *
            >>> m = Grid2D(nx=100, ny=100)
            >>> x, y = m.cellCenters
            >>> v = CellVariable(mesh=m, value=x**2 + y**2 - 10**2)
            >>> s = CellVariable(mesh=m, value=sin(x / 10) * cos(y / 30))
            >>> viewer = MatplotlibSurfactantViewer(distanceVar=v, surfactantVar=s)
            >>> for r in range(1,200):
            ...     v.setValue(x**2 + y**2 - r**2)
            ...     viewer.plot()

            >>> from fipy import *
            >>> dx = 1.
            >>> dy = 1.
            >>> nx = 11
            >>> ny = 11
            >>> Lx = ny * dy
            >>> Ly = nx * dx
            >>> mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)
            >>> # from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
            >>> var = DistanceVariable(mesh = mesh, value = -1)

            >>> x, y = mesh.cellCenters

            >>> var.setValue(1, where=(x - Lx / 2.)**2 + (y - Ly / 2.)**2 < (Lx / 4.)**2)
            >>> var.calcDistanceFunction()
            >>> viewer = MatplotlibSurfactantViewer(var, smooth = 2)
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer

            >>> var = DistanceVariable(mesh = mesh, value = -1)

            >>> var.setValue(1, where=(y > 2. * Ly / 3.) | ((x > Lx / 2.) & (y > Ly / 3.)) | ((y < Ly / 6.) & (x > Lx / 2)))
            >>> var.calcDistanceFunction()
            >>> viewer = MatplotlibSurfactantViewer(var)
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer

            >>> viewer = MatplotlibSurfactantViewer(var, smooth = 2)
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
            >>> del viewer

        :Parameters:

          - `distanceVar`: a `DistanceVariable` object.
          - `levelSetValue`: the value of the contour to be displayed
          - `title`: displayed at the top of the `Viewer` window
          - `animate`: whether to show only the initial condition and the
          - `limits`: a dictionary with possible keys `xmin`, `xmax`,
            `ymin`, `ymax`, `zmin`, `zmax`, `datamin`, `datamax`.  A 1D
            `Viewer` will only use `xmin` and `xmax`, a 2D viewer will also
            use `ymin` and `ymax`, and so on.  All viewers will use
            `datamin` and `datamax`.  Any limit set to a (default) value of
            `None` will autoscale.
            moving top boundary or to show all contours (Default)
        """

        kwlimits.update(limits)
        AbstractMatplotlibViewer.__init__(self, vars=[], title=title, **kwlimits)

        self.distanceVar = distanceVar
        if surfactantVar is None:
            self.surfactantVar = numerix.zeros(len(self.distanceVar), 'd')
        else:
            self.surfactantVar = surfactantVar
        self.smooth = smooth
        self.zoomFactor = zoomFactor

        self.animate = animate
        if animate:
            self._initialCondition = None

        if distanceVar.mesh.dim != 2:
            raise MeshDimensionError('The MatplotlibSurfactantViewer only works for 2D meshes.')

    def _plot(self):
        mesh = self.distanceVar.mesh
        shape = mesh.shape
        X, Y = mesh.cellCenters

        maxX = max(X)

        X = X.reshape(shape, order="FORTRAN")
        Y = Y.reshape(shape, order="FORTRAN")
        Z = self.distanceVar.value.reshape(shape, order="FORTRAN")

        zmin, zmax = self._autoscale(vars=(self.surfactantVar,),
                                     datamin=self._getLimit(('datamin', 'zmin')),
                                     datamax=self._getLimit(('datamax', 'zmax')))

        import pylab
        import matplotlib

        CS = pylab.contour(X, Y, Z, (0.,))
        zc = CS.collections[0]

        verts = numerix.array(zc.get_verts())
        IDs = numerix.array([mesh._getNearestCellID(vert[...,numerix.newaxis]) for vert in verts])

        colors = pylab.cm.jet(( self.surfactantVar[IDs] - zmin) / (zmax - zmin))
        segments = zip(verts[:-1], verts[1:])
        LC = matplotlib.collections.LineCollection(segments, colors=colors)

#         CS.ax.add_collection(LC)

        verts = numerix.array((-verts[...,0], verts[..., 1])).swapaxes(0,1)
        segments = zip(verts[:-1], verts[1:])
        LC = matplotlib.collections.LineCollection(segments, colors=colors)

#         CS.ax.add_collection(LC)

        CS.ax.set_xlim((-maxX, maxX))

#         zc.set_color(pylab.cm.jet(IDs / mesh.numberOfCells))

#         pylab.xlim(xmin=self._getLimit('xmin'),
#                    xmax=self._getLimit('xmax'))
#
#         pylab.ylim(ymin=self._getLimit('ymin'),
#                    ymax=self._getLimit('ymax'))

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
