#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "matplotlibVectorViewer.py"
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

__all__ = ["MatplotlibVectorViewer"]

class MatplotlibVectorViewer(AbstractMatplotlib2DViewer):
    """Displays a vector plot of a 2D rank-1 `CellVariable` or
    `FaceVariable` object using Matplotlib_

    .. _Matplotlib: http://matplotlib.sourceforge.net/

    """

    __doc__ += AbstractMatplotlib2DViewer._test2Dvector(viewer="MatplotlibVectorViewer")
    __doc__ += """

            >>> for sparsity in numerix.arange(5000, 0, -500):
            ...     viewer.quiver(sparsity=sparsity)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()

    """
    __doc__ += AbstractMatplotlib2DViewer._test2DvectorIrregular(viewer="MatplotlibVectorViewer")

    def __init__(self, vars, title=None, scale=None, sparsity=None, log=False, limits={}, axes=None, figaspect='auto', **kwlimits):
        """Creates a `Matplotlib2DViewer`.

        :Parameters:
          vars
            a rank-1 `CellVariable` or `FaceVariable` object.
          title
            displayed at the top of the `Viewer` window
          scale
            if not `None`, scale all arrow lengths by this value
          sparsity
            if not `None`, then this number of arrows will be
            randomly chosen (weighted by the cell volume or face area)
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
        """
        kwlimits.update(limits)
        AbstractMatplotlib2DViewer.__init__(self, vars=vars, title=title, axes=axes, figaspect=figaspect, **kwlimits)

        self.quiver(sparsity=sparsity, scale=scale)
        self.log = log

        self._plot()

    def quiver(self, sparsity=None, scale=None):
        var = self.vars[0]
        mesh = var.mesh

        if isinstance(var, FaceVariable):
            N = mesh.numberOfFaces
            X, Y = mesh.faceCenters

        elif isinstance(var, CellVariable):
            N = mesh.numberOfCells
            X, Y = mesh.cellCenters

        if sparsity is not None and N > sparsity:
            XYrand = numerix.random.random((2, sparsity))
            XYrand = numerix.array([[min(X)],
                                    [min(Y)]]) + XYrand * numerix.array([[max(X) - min(X)],
                                                                         [max(Y) - min(Y)]])
            self.indices = numerix.nearest(numerix.array([X, Y]), XYrand)
        else:
            self.indices = numerix.arange(N)

        X = numerix.take(X, self.indices)
        Y = numerix.take(Y, self.indices)

        U = V = numerix.ones(X.shape, 'l')

        if hasattr(self, "_quiver"):
            self._quiver.remove()

        self._quiver = self.axes.quiver(X, Y, U, V, scale=scale, pivot='middle')

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
            raise MeshDimensionError, "The mesh must be a Mesh2D instance"
        # this viewer can only display one variable
        return [vars[0]]

    def _plot(self):

        var = self.vars[0]
        mesh = var.mesh

        U, V = var.numericValue

        U = numerix.take(U, self.indices)
        V = numerix.take(V, self.indices)

        ang = numerix.arctan2(V, U)
        mag = numerix.sqrt(U**2 + V**2)

        datamin, datamax = self._autoscale(vars=(mag,),
                                           datamin=self._getLimit('datamin'),
                                           datamax=self._getLimit('datamax'))

        mag = numerix.where(mag > datamax, datamax, mag)
        mag = numerix.ma.masked_array(mag, mag < datamin)

        if self.log:
            mag = numerix.log10(mag)
            mag = numerix.ma.masked_array(mag, numerix.isnan(mag))

        U = mag * numerix.cos(ang)
        V = mag * numerix.sin(ang)

        self._quiver.set_UVC(U, V)

        self.axes.set_xlim(xmin=self._getLimit('xmin'),
                           xmax=self._getLimit('xmax'))
        self.axes.set_ylim(ymin=self._getLimit('ymin'),
                           ymax=self._getLimit('ymax'))

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
