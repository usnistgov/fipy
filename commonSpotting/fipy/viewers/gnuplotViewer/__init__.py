__docformat__ = 'restructuredtext'

from gnuplot1DViewer import Gnuplot1DViewer
from gnuplot2DViewer import Gnuplot2DViewer

__all__ = ["GnuplotViewer", "Gnuplot1DViewer", "Gnuplot2DViewer"]

def GnuplotViewer(vars, title=None, limits={}, **kwlimits):
    """Generic function for creating a `GnuplotViewer`. 
    
    The `GnuplotViewer` factory will search the module tree and return an
    instance of the first `GnuplotViewer` it finds of the correct dimension.
    
    :Parameters:
      vars
        a `CellVariable` or tuple of `CellVariable` objects to plot
      title
        displayed at the top of the Viewer window
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

    mesh = vars[0].getMesh()
    
    for var in vars:
        assert mesh is var.getMesh()

    dim = mesh.getDim()
    
    if dim == 1:
        return Gnuplot1DViewer(vars=vars, title=title, **kwlimits)
    elif dim == 2:
        return Gnuplot2DViewer(vars=vars, title=title, **kwlimits)
    else:
        raise IndexError, "Gnuplot can only plot 1D and 2D data"

def make(*args, **kwargs):
    """
    A deprecated synonym for `GnuplotViewer`
    """
    import warnings
    warnings.warn("'GnuplotViewer' should be used instead of 'make'", DeprecationWarning, stacklevel=2)
    return GnuplotViewer(*args, **kwargs)
