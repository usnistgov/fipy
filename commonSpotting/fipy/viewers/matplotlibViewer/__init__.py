__docformat__ = 'restructuredtext'

from matplotlib1DViewer import Matplotlib1DViewer
from matplotlib2DGridViewer import Matplotlib2DGridViewer
from matplotlib2DGridContourViewer import Matplotlib2DGridContourViewer
from matplotlib2DViewer import Matplotlib2DViewer
from matplotlibVectorViewer import MatplotlibVectorViewer

__all__ = ["MatplotlibViewer", "Matplotlib1DViewer", "Matplotlib2DGridViewer", "Matplotlib2DGridContourViewer", "Matplotlib2DViewer", "MatplotlibVectorViewer"]

def MatplotlibViewer(vars, title=None, limits={}, **kwlimits):
    """Generic function for creating a `MatplotlibViewer`. 
    
    The `MatplotlibViewer` factory will search the module tree and return an
    instance of the first `MatplotlibViewer` it finds of the correct dimension
    and rank.
    
    :Parameters:
      vars
        a `CellVariable` or tuple of `CellVariable` objects to plot
      title
        displayed at the top of the `Viewer` window
      limits : dict
        a (deprecated) alternative to limit keyword arguments
      xmin, xmax, ymin, ymax, datamin, datamax
        displayed range of data. A 1D `Viewer` will only use `xmin` and
        `xmax`, a 2D viewer will also use `ymin` and `ymax`. All
        viewers will use `datamin` and `datamax`. Any limit set to a
        (default) value of `None` will autoscale.

    """
    if type(vars) not in [type([]), type(())]:
        vars = [vars]
        
    kwlimits.update(limits)
    
    from fipy.viewers import MeshDimensionError
    
    try:
        return Matplotlib1DViewer(vars=vars, title=title, **kwlimits)
    except MeshDimensionError:
        try:
            from matplotlib2DGridViewer import Matplotlib2DGridViewer
            return Matplotlib2DGridViewer(vars=vars, title=title, **kwlimits)
        except MeshDimensionError:
            try:
                from matplotlib2DViewer import Matplotlib2DViewer
                return Matplotlib2DViewer(vars=vars, title=title, **kwlimits)
            except MeshDimensionError:
                from matplotlibVectorViewer import MatplotlibVectorViewer
                return MatplotlibVectorViewer(vars=vars, title=title, **kwlimits)

def make(*args, **kwargs):
    """
    A deprecated synonym for `MatplotlibViewer`
    """
    import warnings
    warnings.warn("'MatplotlibViewer' should be used instead of 'make'", DeprecationWarning, stacklevel=2)
    return MatplotlibViewer(*args, **kwargs)
