__docformat__ = 'restructuredtext'

from mayaviScalarViewer import MayaviScalarViewer
from mayaviVectorViewer import MayaviVectorViewer

__all__ = ["MayaviViewer", "MayaviScalarViewer", "MayaviVectorViewer"]

def MayaviViewer(vars, title=None, limits={}, **kwlimits):
    """Generic function for creating a `MayaviViewer`. 
    
    The `MayaviViewer` factory will search the module tree and return an
    instance of the first `MayaviViewer` it finds of the correct dimension
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
        return MayaviScalarViewer(vars=vars, title=title, **kwlimits)
    except MeshDimensionError:
        return MayaviVectorViewer(vars=vars, title=title, **kwlimits)

def make(*args, **kwargs):
    """
    A deprecated synonym for `MayaviViewer`
    """
    import warnings
    warnings.warn("'MayaviViewer' should be used instead of 'make'", DeprecationWarning, stacklevel=2)
    return MatplotlibViewer(*args, **kwargs)
