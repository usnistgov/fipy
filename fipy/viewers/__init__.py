__docformat__ = 'restructuredtext'

import os
import sys
import string
import glob
import imp

try:
    from gistViewer import *
except:
    pass

try:
    from gnuplotViewer import *
except:
    pass

try:
    from matplotlibViewer import *
except:
    pass

try:
    from mayaviViewer import *
except:
    pass

from multiViewer import MultiViewer
from tsvViewer import TSVViewer
from vtkViewer import VTKViewer


# what about vector variables?

class MeshDimensionError(IndexError):
    pass
    
from viewer import _Viewer
class DummyViewer(_Viewer):
    def plot(self, filename=None):
        pass

def Viewer(vars, title=None, limits={}, FIPY_VIEWER=None, **kwlimits):
    r"""Generic function for creating a `Viewer`. 
    
    The `Viewer` factory will search the module tree and return an instance of
    the first `Viewer` it finds that supports the dimensions of `vars`. Setting
    the '`FIPY_VIEWER`' environment variable to either '`gist`', '`gnuplot`',
    '`matplotlib`', '`tsv`', or '`vtk`' will specify the viewer.
       
    The `kwlimits` or `limits` parameters can be used to constrain the view. For example::
            
        Viewer(vars=some1Dvar, xmin=0.5, xmax=None, datamax=3)
               
    or::
        
        Viewer(vars=some1Dvar, 
               limits={'xmin': 0.5, 'xmax': None, 'datamax': 3})
        
    will return a viewer that displays a line plot from an `x` value
    of 0.5 up to the largest `x` value in the dataset. The data values
    will be truncated at an upper value of 3, but will have no lower
    limit.

    :Parameters:
      vars
        a `CellVariable` or tuple of `CellVariable` objects to plot
      title
        displayed at the top of the `Viewer` window
      limits : dict
        a (deprecated) alternative to limit keyword arguments
      FIPY_VIEWER
        a specific viewer to attempt (possibly multiple times for multiple variables)
      xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax
        displayed range of data. A 1D `Viewer` will only use `xmin` and
        `xmax`, a 2D viewer will also use `ymin` and `ymax`, and so on. All
        viewers will use `datamin` and `datamax`. Any limit set to a
        (default) value of `None` will autoscale.
      
    """
    
    if type(vars) not in [type([]), type(())]:
        vars = [vars]
    vars = list(vars)
    
    if FIPY_VIEWER is None and os.environ.has_key('FIPY_VIEWER'):
        FIPY_VIEWER = os.environ['FIPY_VIEWER']

    if FIPY_VIEWER == "dummy":
        return DummyViewer(vars=vars)

    errors = []

    attempts = []
    viewers = []
    
    import pkg_resources
    for ep in pkg_resources.iter_entry_points(group='fipy.viewers', 
                                              name=FIPY_VIEWER):
                                                  
        attempts.append(ep.name)
        
        try:
            ViewerClass = ep.load()
            
            while len(vars) > 0:
                viewer = ViewerClass(vars=vars, title=title, limits=limits, **kwlimits)
                
                for var in viewer.getVars():
                    vars.remove(var)
                
                viewers.append(viewer)
                
            break
        except Exception, s:
            errors.append("%s: %s" % (ep.name, s))

    if len(attempts) == 0:
        if FIPY_VIEWER is not None:
            raise ImportError, "`%s` viewer not found" % FIPY_VIEWER
        else:
            raise ImportError, "No viewers found"
    
    if len(vars) > 0:
        raise ImportError, "Failed to import a viewer: %s" % str(errors)
            
    if len(viewers) > 1:
        return MultiViewer(viewers = viewers)
    else:
        return viewers[0]
        
def make(*args, **kwargs):
    """
    A deprecated synonym for `Viewer`
    """
    import warnings
    warnings.warn("'Viewer' should be used instead of 'make'", DeprecationWarning, stacklevel=2)
    return Viewer(*args, **kwargs)
