from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.viewers.matplotlibViewer.matplotlib1DViewer import *
from fipy.viewers.matplotlibViewer.matplotlib2DGridViewer import *
from fipy.viewers.matplotlibViewer.matplotlib2DContourViewer import *
from fipy.viewers.matplotlibViewer.matplotlib2DGridContourViewer import *
from fipy.viewers.matplotlibViewer.matplotlib2DViewer import *
from fipy.viewers.matplotlibViewer.matplotlibVectorViewer import *
from fipy.viewers.matplotlibViewer.matplotlibStreamViewer import *

def MatplotlibViewer(vars, title=None, limits={}, cmap=None, colorbar='vertical', axes=None, **kwlimits):
    """Generic function for creating a `MatplotlibViewer`.

    The `MatplotlibViewer` factory will search the module tree and return an
    instance of the first `MatplotlibViewer` it finds of the correct dimension
    and rank.

    It is possible to view different
    :class:`~fipy.variables.variable.Variable`\s against different Matplotlib_
    `Axes`

    >>> from matplotlib import pyplot as plt
    >>> from fipy import *

    >>> plt.ion()
    >>> fig = plt.figure()

    >>> ax1 = plt.subplot((221))
    >>> ax2 = plt.subplot((223))
    >>> ax3 = plt.subplot((224))

    >>> k = Variable(name="k", value=0.)

    >>> mesh1 = Grid1D(nx=100)
    >>> x, = mesh1.cellCenters
    >>> xVar = CellVariable(mesh=mesh1, name="x", value=x)
    >>> viewer1 = MatplotlibViewer(vars=(numerix.sin(0.1 * k * xVar), numerix.cos(0.1 * k * xVar / numerix.pi)),
    ...                            limits={'xmin': 10, 'xmax': 90},
    ...                            datamin=-0.9, datamax=2.0,
    ...                            title="Grid1D test",
    ...                            axes=ax1,
    ...                            legend=None)

    >>> mesh2 = Grid2D(nx=50, ny=100, dx=0.1, dy=0.01)
    >>> x, y = mesh2.cellCenters
    >>> xyVar = CellVariable(mesh=mesh2, name="x y", value=x * y)
    >>> viewer2 = MatplotlibViewer(vars=numerix.sin(k * xyVar),
    ...                            limits={'ymin': 0.1, 'ymax': 0.9},
    ...                            datamin=-0.9, datamax=2.0,
    ...                            title="Grid2D test",
    ...                            axes=ax2,
    ...                            colorbar=None)

    >>> mesh3 = (Grid2D(nx=5, ny=10, dx=0.1, dy=0.1)
    ...          + (Tri2D(nx=5, ny=5, dx=0.1, dy=0.1)
    ...             + ((0.5,), (0.2,))))
    >>> x, y = mesh3.cellCenters
    >>> xyVar = CellVariable(mesh=mesh3, name="x y", value=x * y)
    >>> viewer3 = MatplotlibViewer(vars=numerix.sin(k * xyVar),
    ...                            limits={'ymin': 0.1, 'ymax': 0.9},
    ...                            datamin=-0.9, datamax=2.0,
    ...                            title="Irregular 2D test",
    ...                            axes=ax3,
    ...                            cmap = plt.cm.OrRd)

    >>> viewer = MultiViewer(viewers=(viewer1, viewer2, viewer3))
    >>> from builtins import range
    >>> for kval in range(10):
    ...     k.setValue(kval)
    ...     viewer.plot()

    >>> viewer._promptForOpinion()

    .. _Matplotlib: http://matplotlib.sourceforge.net/

    Parameters
    ----------
    vars : ~fipy.variables.cellVariable.CellVariable or list
        the `Variable` objects to display.
    title : str, optional
        displayed at the top of the `Viewer` window
    limits : dict
        a (deprecated) alternative to limit keyword arguments
    xmin, xmax, ymin, ymax, datamin, datamax : float, optional
        displayed range of data. A 1D `Viewer` will only use *xmin* and
        *xmax*, a 2D viewer will also use *ymin* and *ymax*. All
        viewers will use *datamin* and *datamax*. Any limit set to a
        (default) value of `None` will autoscale.
    cmap : ~matplotlib.colors.Colormap, optional
        the :class:`~matplotlib.colors.Colormap`.
        Defaults to `matplotlib.cm.jet`
    colorbar : bool, optional
        plot a color bar in specified orientation if not `None`
    axes : ~matplotlib.axes.Axes, optional
        if not `None`, `vars` will be plotted into this Matplotlib `Axes` object

    """
    if type(vars) not in [type([]), type(())]:
        vars = [vars]

    kwlimits.update(limits)

    from fipy.viewers import MeshDimensionError
    from fipy.viewers.matplotlibViewer import (Matplotlib1DViewer,
                                               Matplotlib2DGridViewer,
                                               Matplotlib2DViewer,
                                               MatplotlibVectorViewer)

    try:
        return Matplotlib1DViewer(vars=vars, title=title, axes=axes, **kwlimits)
    except MeshDimensionError:
        try:
            from fipy.viewers.matplotlibViewer.matplotlib2DGridViewer import Matplotlib2DGridViewer
            return Matplotlib2DGridViewer(vars=vars, title=title, cmap=cmap, colorbar=colorbar, axes=axes, **kwlimits)
        except MeshDimensionError:
            try:
                from fipy.viewers.matplotlibViewer.matplotlib2DViewer import Matplotlib2DViewer
                return Matplotlib2DViewer(vars=vars, title=title, cmap=cmap, colorbar=colorbar, axes=axes, **kwlimits)
            except MeshDimensionError:
                from fipy.viewers.matplotlibViewer.matplotlibVectorViewer import MatplotlibVectorViewer
                return MatplotlibVectorViewer(vars=vars, title=title, axes=axes, **kwlimits)
