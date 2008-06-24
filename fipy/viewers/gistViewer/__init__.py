__docformat__ = 'restructuredtext'

from gist1DViewer import Gist1DViewer
from gist2DViewer import Gist2DViewer
from gistVectorViewer import GistVectorViewer

__all__ = ["Gist1DViewer", "Gist2DViewer", "GistVectorViewer"]

def GistViewer(vars, title = None, limits = None):
    r"""
    Generic function for creating a `GistViewer`. The `make` function
    will search the module tree and return an instance of the first
    `GistViewer` it finds of the correct dimension.
        
    :Parameters:

      - `vars`: a `CellVariable` or tuple of `CellVariable` objects to plot
      - `limits`: a dictionary with possible keys `'xmin'`, `'xmax'`,
        `'ymin'`, `'ymax'`, `'zmin'`, `'zmax'`, `'datamin'`, `'datamax'`.
        A 1D `Viewer` will only use `'xmin'` and `'xmax'`, a 2D viewer
        will also use `'ymin'` and `'ymax'`, and so on.
        All viewers will use `'datamin'` and `'datamax'`.
        Any limit set to a (default) value of `None` will autoscale.
      - `title`: displayed at the top of the `Viewer` window
    """
    if type(vars) not in [type([]), type(())]:
        vars = [vars]
        
    from fipy.viewers import MeshDimensionError
    
    try:
        return Gist1DViewer(vars = vars, title = title, limits = limits)
    except MeshDimensionError:
        try:
            return Gist2DViewer(vars = vars, title = title, limits = limits)
        except MeshDimensionError:
            return GistVectorViewer(vars = vars, title = title)
            
def make(*args, **kwargs):
    import warnings
    warnings.warn("'GistViewer' should be used instead of 'make'", DeprecationWarning, stacklevel=2)
    return GistViewer(*args, **kwargs)
