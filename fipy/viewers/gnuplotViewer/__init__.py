from gnuplotViewer import GnuplotViewer

def make(vars, title = None, limits = None):
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
