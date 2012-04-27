__docformat__ = 'restructuredtext'

from fipy.viewers.vtkViewer.vtkCellViewer import VTKCellViewer
from fipy.viewers.vtkViewer.vtkFaceViewer import VTKFaceViewer

__all__ = ["VTKViewer"]
__all__.extend(vtkCellViewer.__all__)
__all__.extend(vtkFaceViewer.__all__)

def VTKViewer(vars, title=None, limits={}, **kwlimits):
    """Generic function for creating a `VTKViewer`. 
    
    The `VTKViewer` factory will search the module tree and return an
    instance of the first `VTKViewer` it finds of the correct dimension
    and rank.
    
    :Parameters:
      vars
        a `_MeshVariable` or tuple of `_MeshVariable` objects to plot
      title
        displayed at the top of the `Viewer` window
      limits : dict
        a (deprecated) alternative to limit keyword arguments
      xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax
        displayed range of data. A 1D `Viewer` will only use `xmin` and
        `xmax`, a 2D viewer will also use `ymin` and `ymax`, and so on. All
        viewers will use `datamin` and `datamax`. Any limit set to a
        (default) value of `None` will autoscale.

    """
    if type(vars) not in [type([]), type(())]:
        vars = [vars]
        
    kwlimits.update(limits)
    
    try:
        return VTKCellViewer(vars=vars, title=title, **kwlimits)
    except TypeError:
        return VTKFaceViewer(vars=vars, title=title, **kwlimits)
