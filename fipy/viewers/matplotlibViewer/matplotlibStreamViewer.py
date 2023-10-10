from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.variables.faceVariable import FaceVariable
from fipy.variables.cellVariable import CellVariable

from fipy.viewers.matplotlibViewer.abstractMatplotlib2DViewer import AbstractMatplotlib2DViewer

__all__ = ["MatplotlibStreamViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

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

    def __init__(self, vars, title=None, log=False, limits={}, axes=None, figaspect='auto',
                 density=1, linewidth=None, color=None, cmap=None, norm=None, arrowsize=1,
                 arrowstyle='-|>', minlength=0.1,
                 **kwlimits):
        """Creates a `MatplotlibStreamViewer`.

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable or ~fipy.variables.faceVariable.FaceVariable
            rank-1 `Variable` to display
        title : str, optional
            displayed at the top of the `Viewer` window
        log : bool, optional
            if `True`, arrow length goes at the base-10 logarithm of the magnitude
        limits : dict, optional
            a (deprecated) alternative to limit keyword arguments
        xmin, xmax, ymin, ymax, datamin, datamax : float
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
        axes : ~matplotlib.axes.Axes, optional
            if not `None`, `vars` will be plotted into this Matplotlib `Axes` object
        figaspect : float, optional
            desired aspect ratio of figure. If arg is a number, use that aspect
            ratio. If arg is `auto`, the aspect ratio will be determined from
            the Variable's mesh.
        density : float or tuple of float, optional
            Controls the closeness of streamlines.  When ``density = 1``,
            the domain is divided into a 30x30 grid.  *density* linearly
            scales this grid.  Each cell in the grid can have, at most, one
            traversing streamline.  For different densities in each
            direction, use a tuple (density_x, density_y).
        linewidth : array_like or ~fipy.variables.cellVariable.CellVariable or ~fipy.variables.faceVariable.FaceVariable, optional
            The width of the stream lines.  With a rank-0 `CellVariable` or
            `FaceVariable` the line width can be varied across the grid.
            The MeshVariable must have the same type and be defined on
            the same `Mesh` as *vars*.
        color : str or ~fipy.variables.cellVariable.CellVariable or ~fipy.variables.faceVariable.FaceVariable, optional
            The streamline color as a matplotlib color code or a field of
            numbers.  If given a rank-0 `CellVariable` or `FaceVariable`,
            its values are converted to colors using *cmap* and *norm*.
            The MeshVariable must have the same type and be defined on the
            same `Mesh` as *vars*.
        cmap : ~matplotlib.colors.Colormap, optional
            Colormap used to plot streamlines and arrows.  This is only
            used if *color* is a MeshVariable.
        norm : ~matplotlib.colors.Normalize, optional
            Normalize object used to scale luminance data to 0, 1.  If
            ``None``, stretch (min, max) to (0, 1).  Only necessary when
            *color* is a MeshVariable.
        arrowsize : float, optional
            Scaling factor for the arrow size.
        arrowstyle : str, optional
            Arrow style specification.
            See `~matplotlib.patches.FancyArrowPatch`.
        minlength : float, optional
              Minimum length of streamline in axes coordinates.

        """
        kwlimits.update(limits)
        AbstractMatplotlib2DViewer.__init__(self, vars=vars, title=title, axes=axes, figaspect=figaspect, **kwlimits)

        self.log = log
        self.kwargs = dict(density=density, linewidth=linewidth, color=color,
                           cmap=cmap, norm=norm, arrowsize=arrowsize,
                           arrowstyle=arrowstyle, minlength=minlength)

        self._stream = None

        self._plot()

    @property
    def kwargs(self):
        """keyword arguments to pass to :meth:`~matplotlib.axes.Axes.streamplot`."""
        return self._kwargs

    @kwargs.setter
    def kwargs(self, value):
        self._kwargs = value

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

        lw = self.kwargs["linewidth"]
        if isinstance(lw, (FaceVariable, CellVariable)):
            lw = griddata(C.value.T, lw.value,
                          (grid_x, grid_y), method='cubic')

        color = self.kwargs["color"]
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

        kwargs = self.kwargs.copy()
        kwargs["linewidth"] = lw
        kwargs["color"] = color

        self.axes.cla()
        self._stream = self.axes.streamplot(X, Y, U, V, **kwargs)

        self.axes.set_xlim(xmin=self._getLimit('xmin'),
                           xmax=self._getLimit('xmax'))
        self.axes.set_ylim(ymin=self._getLimit('ymin'),
                           ymax=self._getLimit('ymax'))

    @classmethod
    def _doctest_body(cls):
        return (cls._test2Dvector()
                + cls._test2DvectorIrregular())

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
