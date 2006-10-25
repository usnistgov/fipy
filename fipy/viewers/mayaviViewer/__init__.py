__docformat__ = 'restructuredtext'

from mayaviViewer import MayaviViewer
from mayaviDistanceViewer import MayaviDistanceViewer
from mayaviSurfactantViewer import MayaviSurfactantViewer

__all__ = ["MayaviViewer", "MayaviDistanceViewer", "MayaviSurfactantViewer"]

def make(vars, title = None, limits = None):
    r"""
    Generic function for creating a `MayaviViewer`..
        
    :Parameters:

      - `vars`: a `CellVariable` or tuple of `CellVariable` objects to plot
      - `limits`: a dictionary with possible keys `'xmin'`, `'xmax'`,
        `'ymin'`, `'ymax'`, `'zmin'`, `'zmax'`, `'datamin'`, `'datamax'`.
        A 1D Viewer will only use `'xmin'` and `'xmax'`, a 2D viewer
        will also use `'ymin'` and `'ymax'`, and so on.
        All viewers will use `'datamin'` and `'datamax'`.
        Any limit set to a (default) value of `None` will autoscale.
      - `title`: displayed at the top of the Viewer window
    """
    if type(vars) not in [type([]), type(())]:
        vars = [vars]

    for var in vars:
        if var.getMesh().getDim() == 1:
            raise IndexError, "Mayavi does not display useful 1D data"

    return MayaviViewer(vars, title = title, limits = limits) 
        

