__docformat__ = 'restructuredtext'

__all__ = []

try:
    from fipy.viewers.matplotlibViewer import *
    __all__.extend(matplotlibViewer.__all__)
except:
    pass

try:
    from fipy.viewers.mayaviViewer import *
    __all__.extend(mayaviViewer.__all__)
except:
    pass

from fipy.viewers.multiViewer import *
from fipy.viewers.tsvViewer import *
from fipy.viewers.vtkViewer import *

__all__.extend(multiViewer.__all__)
__all__.extend(tsvViewer.__all__)
__all__.extend(vtkViewer.__all__)

# what about vector variables?

class MeshDimensionError(IndexError):
    pass

from fipy.viewers.viewer import AbstractViewer
class DummyViewer(AbstractViewer):
    def plot(self, filename=None):
        pass

def Viewer(vars, title=None, limits={}, FIPY_VIEWER=None, **kwlimits):
    r"""Generic function for creating a `Viewer`.

    The `Viewer` factory will search the module tree and return an instance of
    the first `Viewer` it finds that supports the dimensions of `vars`. Setting
    the '`FIPY_VIEWER`' environment variable to either '`matplotlib`', '`mayavi`',
    '`tsv`', or '`vtk`' will specify the viewer.

    The `kwlimits` or `limits` parameters can be used to constrain the view. For example::

        Viewer(vars=some1Dvar, xmin=0.5, xmax=None, datamax=3)

    or::

        Viewer(vars=some1Dvar,
               limits={'xmin': 0.5, 'xmax': None, 'datamax': 3})

    will return a viewer that displays a line plot from an `x` value
    of 0.5 up to the largest `x` value in the dataset. The data values
    will be truncated at an upper value of 3, but will have no lower
    limit.

    :Parameters:
      vars
        a `CellVariable` or tuple of `CellVariable` objects to plot
      title
        displayed at the top of the `Viewer` window
      limits : dict
        a (deprecated) alternative to limit keyword arguments
      FIPY_VIEWER
        a specific viewer to attempt (possibly multiple times for multiple variables)
      xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax
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

    enpts = []
    import pkg_resources
    for ep in pkg_resources.iter_entry_points(group='fipy.viewers',
                                              name=FIPY_VIEWER):
        enpts.append((ep.name,ep))

    for name, ep in sorted(enpts):

        attempts.append(name)

        try:
            ViewerClass = ep.load()

            while len(vars) > 0:
                viewer = ViewerClass(vars=vars, title=title, limits=limits, **kwlimits)

                for var in viewer.vars:
                    vars.remove(var)

                viewers.append(viewer)

            break
        except Exception as s:
            errors.append("%s: %s" % (name, s))

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

__all__.extend(["MeshDimensionError", "DummyViewer", "Viewer"])
