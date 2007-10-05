__docformat__ = 'restructuredtext'

from matplotlib1DViewer import Matplotlib1DViewer
from matplotlib2DGridViewer import Matplotlib2DGridViewer
from matplotlib2DGridContourViewer import Matplotlib2DGridContourViewer
from matplotlib2DViewer import Matplotlib2DViewer
from matplotlibVectorViewer import MatplotlibVectorViewer

__all__ = ["Matplotlib1DViewer", "Matplotlib2DGridViewer", "Matplotlib2DGridContourViewer", "Matplotlib2DViewer", "MatplotlibVectorViewer"]

def make(vars, title = None, limits = None):
    """
    Generic function for creating a `MatplotlibViewer`. The `make` function
    will search the module tree and return an instance of the first
    `MatplotlibViewer` it finds of the correct dimension.
    
    :Parameters:

      - `vars`: a `CellVariable`, or `FaceVariable` object 
        or sequence to plot
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
        
    from fipy.viewers import MeshDimensionError
    
    try:
        return Matplotlib1DViewer(vars = vars, title = title, limits = limits)
    except MeshDimensionError:
##        try:
##            from matplotlib2DGridViewer import Matplotlib2DGridViewer
##            return Matplotlib2DGridViewer(vars=vars, title=title, limits=limits)
##        except MeshDimensionError:
        try:
            from matplotlib2DViewer import Matplotlib2DViewer
            return Matplotlib2DViewer(vars = vars, title = title, limits = limits)
        except MeshDimensionError:
            from matplotlibVectorViewer import MatplotlibVectorViewer
            return MatplotlibVectorViewer(vars = vars, title = title)
