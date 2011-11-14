__docformat__ = 'restructuredtext'

from fipy.viewers.gnuplotViewer.gnuplot1DViewer import *
from fipy.viewers.gnuplotViewer.gnuplot2DViewer import *

__all__ = ["GnuplotViewer"]
__all__.extend(gnuplot1DViewer.__all__)
__all__.extend(gnuplot2DViewer.__all__)

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

    mesh = vars[0].mesh
    
    for var in vars:
        assert mesh is var.mesh

    dim = mesh.dim
    
    if dim == 1:
        return Gnuplot1DViewer(vars=vars, title=title, **kwlimits)
    elif dim == 2:
        return Gnuplot2DViewer(vars=vars, title=title, **kwlimits)
    else:
        raise IndexError, "Gnuplot can only plot 1D and 2D data"
