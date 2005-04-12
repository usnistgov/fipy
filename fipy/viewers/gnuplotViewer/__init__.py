__docformat__ = 'restructuredtext'

from gnuplotViewer import _GnuplotViewer

def make(vars, title = None, limits = None):
    """
    Generic function for creating a `GnuplotViewer`. The `make` function
    will search the module tree and return an instance of the first
    `GnuplotViewer` it finds of the correct dimension.
    Usage:

    ::

        viewer = make(vars)
        viewer.plot()
    
    :Parameters:

      - `vars`: a `CellVariable` or tuple of `CellVariable` objects to plot
      - `limits`: a dictionary with possible keys `xmin`, `xmax`,
        `ymin`, `ymax`, `zmin`, `zmax`, `datamin`, `datamax`.
        A 1D Viewer will only use `xmin` and `xmax`, a 2D viewer
        will also use `ymin` and `ymax`, and so on.
        All viewers will use `datamin` and `datamax`.
        Any limit set to a (default) value of `None` will autoscale.
      - `title`: displayed at the top of the Viewer window

    """
    if type(vars) not in [type([]), type(())]:
        vars = [vars]
        
    mesh = vars[0].getMesh()
    
    for var in vars:
        assert mesh is var.getMesh()

    dim = mesh.getDim()
    
    if dim == 1:
        from gnuplot1DViewer import Gnuplot1DViewer
        return Gnuplot1DViewer(vars = vars, title = title, limits = limits)
    elif dim == 2:
        from gnuplot2DViewer import Gnuplot2DViewer
        return Gnuplot2DViewer(vars = vars, title = title, limits = limits)
    else:
        raise IndexError, "Gist can only plot 1D and 2D data"
