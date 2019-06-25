from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.viewers.vtkViewer.vtkCellViewer import VTKCellViewer
from fipy.viewers.vtkViewer.vtkFaceViewer import VTKFaceViewer

__all__ = ["VTKViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]
__all__.extend(vtkCellViewer.__all__)
__all__.extend(vtkFaceViewer.__all__)

def VTKViewer(vars, title=None, limits={}, **kwlimits):
    """Generic function for creating a `VTKViewer`.

    The `VTKViewer` factory will search the module tree and return an
    instance of the first `VTKViewer` it finds of the correct dimension
    and rank.

    Parameters
    ----------
    vars : ~fipy.variables.cellVariable.CellVariable or ~fipy.variables.faceVariable.FaceVariable or list
        the `MeshVariable` objects to display.
    title : str, optional
        displayed at the top of the `Viewer` window
    limits : dict, optional
        a (deprecated) alternative to limit keyword arguments
    float xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax : float, optional
        displayed range of data. Any limit set to
        a (default) value of `None` will autoscale.
    """
    if type(vars) not in [type([]), type(())]:
        vars = [vars]

    kwlimits.update(limits)

    try:
        return VTKCellViewer(vars=vars, title=title, **kwlimits)
    except TypeError:
        return VTKFaceViewer(vars=vars, title=title, **kwlimits)
