from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from future.builtins import super

from fipy.tools import numerix
from fipy.variables.faceVariable import FaceVariable
from fipy.variables.cellVariable import CellVariable

from fipy.viewers.matplotlibViewer.abstractMatplotlib2DViewer import AbstractMatplotlib2DViewer

__all__ = ["MatplotlibVectorViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class MatplotlibVectorViewer(AbstractMatplotlib2DViewer):
    """Displays a vector plot of a 2D rank-1 :class:`~fipy.variables.meshVariable.MeshVariable` using Matplotlib_

    .. _Matplotlib: http://matplotlib.sourceforge.net/

    """

    def __init__(self, vars, title=None, scale=None, sparsity=None, log=False, limits={}, axes=None, figaspect='auto', **kwlimits):
        """Creates a `Matplotlib2DViewer`.

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable or ~fipy.variables.faceVariable.FaceVariable
            rank-1 `Variable` to display
        title : str, optional
            displayed at the top of the `Viewer` window
        scale : float, optional
            if not `None`, scale all arrow lengths by this value
        sparsity : int, optional
            if not `None`, then this number of arrows will be
            randomly chosen (weighted by the cell volume or face area)
        log : bool, optional
            if `True`, arrow length goes at the base-10 logarithm of the magnitude
        limits : dict, optional
            a (deprecated) alternative to limit keyword arguments
        xmin, xmax, ymin, ymax, datamin, datamax : float, optional
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
        axes : ~matplotlib.axes.Axes, optional
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
        figaspect : float, optional
            desired aspect ratio of figure. If arg is a number, use that aspect
            ratio. If arg is `auto`, the aspect ratio will be determined from
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
            raise MeshDimensionError("The mesh must be a Mesh2D instance")
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

    @classmethod
    def _doctest_body(cls):
        return (cls._test2Dvector()
                + cls._test2DvectorIrregular())

    @classmethod
    def _doctest_extra(cls):
        return ("""
            >>> for sparsity in numerix.arange(5000, 0, -500):
            ...     viewer.quiver(sparsity=sparsity)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()
        """ + super()._doctest_extra())

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
