"""Tools for displaying the values of :class:`~fipy.variables.variable.Variable` objects
"""

from __future__ import unicode_literals
from builtins import str
__docformat__ = 'restructuredtext'

try:
    from fipy.viewers.matplotlibViewer import *
except:
    pass

try:
    from fipy.viewers.mayaviViewer import *
except:
    pass

from fipy.viewers.multiViewer import *
from fipy.viewers.tsvViewer import *
from fipy.viewers.vtkViewer import *

# what about vector variables?

class MeshDimensionError(IndexError):
    pass

from fipy.viewers.viewer import AbstractViewer
class DummyViewer(AbstractViewer):
    """Substitute viewer that doesn't do anything"""
    def plot(self, filename=None):
        pass

def Viewer(vars, title=None, limits={}, FIPY_VIEWER=None, **kwlimits):
    r"""Generic function for creating a `Viewer`.

    The `Viewer` factory will search the module tree and return an instance of
    the first `Viewer` it finds that supports the dimensions of `vars`. Setting
    the `FIPY_VIEWER` environment variable to either `matplotlib`, `mayavi`,
    `tsv`, or `vtk` will specify the viewer.

    The `kwlimits` or `limits` parameters can be used to constrain the view. For example::

        Viewer(vars=some1Dvar, xmin=0.5, xmax=None, datamax=3)

    or::

        Viewer(vars=some1Dvar,
               limits={'xmin': 0.5, 'xmax': None, 'datamax': 3})

    will return a viewer that displays a line plot from an `x` value
    of 0.5 up to the largest `x` value in the dataset. The data values
    will be truncated at an upper value of 3, but will have no lower
    limit.

    Parameters
    ----------
    vars : ~fipy.variables.cellVariable.CellVariable or list
        the `Variable` objects to display.
    title : str, optional
        displayed at the top of the `Viewer` window
    limits : dict
        a (deprecated) alternative to limit keyword arguments
    FIPY_VIEWER
        a specific viewer to attempt (possibly multiple times for multiple variables)
    xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax : float, optional
        displayed range of data. A 1D `Viewer` will only use `xmin` and
        `xmax`, a 2D viewer will also use `ymin` and `ymax`, and so on. All
        viewers will use `datamin` and `datamax`. Any limit set to a
        (default) value of `None` will autoscale.
    """
    import os

    if type(vars) not in [type([]), type(())]:
        vars = [vars]
    vars = list(vars)

    if FIPY_VIEWER is None and 'FIPY_VIEWER' in os.environ:
        FIPY_VIEWER = os.environ['FIPY_VIEWER']

    if FIPY_VIEWER == "dummy":
        return DummyViewer(vars=vars)

    errors = []

    attempts = []
    viewers = []

    emptyvars = [var for var in vars if var.mesh.numberOfCells == 0]
    vars = [var for var in vars if var.mesh.numberOfCells > 0]

    if len(emptyvars):
        viewers.append(DummyViewer(vars=emptyvars))

    try:
        # pkg_resources is deprecated,
        # but importlib.metadata doesn't exist until Python 3.8
        from importlib.metadata import entry_points

        enpts = entry_points()

        if hasattr(enpts, "select"):
            # .select() not introduced until
            # importlib_metadata 3.6 and Python 3.10

            if FIPY_VIEWER is None:
                # unlike pkg_resources.iter_entry_points,
                # importlib.metadata.entry_points doesn't return anything
                # if name=NONE
                enpts = enpts.select(group='fipy.viewers')
            else:
                enpts = enpts.select(group='fipy.viewers', name=FIPY_VIEWER)
        else:
            enpts = enpts.get("fipy.viewers", ())

            if FIPY_VIEWER is not None:
                enpts = (mod for mod in enpts if mod.name == FIPY_VIEWER)

        enpts = sorted(enpts)
    except ImportError:
        from pkg_resources import iter_entry_points

        enpts = iter_entry_points(group='fipy.viewers', name=FIPY_VIEWER)

        # pkg_resources.EntryPoint objects aren't sortable
        enpts = sorted(enpts, key=lambda ep: ep.name)

    for ep in enpts:

        attempts.append(ep.name)

        try:
            ViewerClass = ep.load()

            while len(vars) > 0:
                viewer = ViewerClass(vars=vars, title=title, limits=limits, **kwlimits)

                for var in viewer.vars:
                    vars.remove(var)

                viewers.append(viewer)

            break
        except Exception as s:
            errors.append("%s: %s" % (ep.name, s))

    if len(attempts) == 0:
        if FIPY_VIEWER is not None:
            raise ImportError("`%s` viewer not found" % FIPY_VIEWER)
        else:
            raise ImportError("No viewers found. Run `python setup.py egg_info` or similar.")

    if len(vars) > 0:
        raise ImportError("Failed to import a viewer: %s" % str(errors))

    if len(viewers) > 1:
        return MultiViewer(viewers = viewers)
    else:
        return viewers[0]
