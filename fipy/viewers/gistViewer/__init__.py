def make(vars, title = None, limits = None):
    if type(vars) not in [type([]), type(())]:
        vars = [vars]
        
    mesh = vars[0].getMesh()
    
    for var in vars:
        assert mesh is var.getMesh()

    dim = mesh.getDim()
    
    if dim == 1:
        from gist1DViewer import Gist1DViewer
        return Gist1DViewer(vars = vars, title = title, limits = limits)
    elif dim == 2:
        from gist2DViewer import Gist2DViewer
        return Gist2DViewer(vars = vars, title = title, limits = limits)
    elif dim == 3:
        raise IndexError, "Gist can only plot 1D and 2D data"
        

