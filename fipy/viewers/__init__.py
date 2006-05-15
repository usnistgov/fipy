__docformat__ = 'restructuredtext'

import os
import sys
import string
import glob
import imp

# what about vector variables?

class MeshDimensionError(IndexError):
    pass

def make(vars, title = None, limits = None):
    r"""
    
    Generic function for creating a `Viewer`. The `make` function will
    search the module tree and return an instance of the first
    `Viewer` it finds that supports the dimensions of `vars`. Setting the
    '`FIPY_VIEWER`' environment variable to either '`gist`',
    '`gnuplot`', '`matplotlib`', or '`tsv`' will specify the viewer.
       
    The `limits` parameter can be used to constrain the view. For example::
            
        fipy.viewers.make(vars = some1Dvar, 
                          limits = {'xmin': 0.5, 'xmax': None, 'datamax': 3})
        
    will return a viewer that displays a line plot from an `x` value
    of 0.5 up to the largest `x` value in the dataset. The data values
    will be truncated at an upper value of 3, but will have no lower
    limit.

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
    vars = list(vars)
    
    if os.environ.has_key('FIPY_VIEWER'):
        viewerClassNames = [os.environ['FIPY_VIEWER'] + 'Viewer']
    else:
        viewerPaths = []
    
        for suffix, mode, moduleType in [("", "", imp.PKG_DIRECTORY)] + imp.get_suffixes():
            if moduleType is not imp.PY_COMPILED:
                for path in __path__:
                    viewerPaths += glob.glob(os.path.join(path, "*Viewer%s" % suffix))
                
        viewerClassNames = []
        for viewerPath in viewerPaths:
            path, f = os.path.split(viewerPath)
            className, ext = os.path.splitext(f)
            viewerClassNames.append(className)


    errors = []

    viewers = []
    for className in viewerClassNames:
        try:
            className = string.lower(className[0]) + className[1:]
            viewerModule = imp.load_module(className, *imp.find_module(className, __path__))
            
            while len(vars) > 0:
                viewer = viewerModule.make(vars = vars, title = title, limits = limits)
                
                for var in viewer.getVars():
                    vars.remove(var)
                
                viewers.append(viewer)
            
            break
        except Exception, s:
            errors.append("%s: %s" % (className, s))
        
    if len(vars) > 0:
        raise ImportError, "Failed to import a viewer: %s" % str(errors)        
            
    if len(viewers) > 1:
        from fipy.viewers.multiViewer import _MultiViewer
        return _MultiViewer(viewers = viewers)
    else:
        return viewers[0]
        

