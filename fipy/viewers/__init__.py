__docformat__ = 'restructuredtext'

import os
import sys
import string
import glob
import imp

# what about vector variables?

def make(vars, title = None, limits = None):
    r"""
    
    Generic function for creating a `Viewer`. The `make` function will
    search the module tree and return an instance of the first
    `Viewer` it finds of the correct dimension. Setting the
    `'FIPY_VIEWER'` environment variable to either `'gist'`,
    `'gnuplot'` or `'tsv'` will specify the viewer.

    ::

       viewer = make(vars)
       viewer.plot()

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
            path, file = os.path.split(viewerPath)
            className, ext = os.path.splitext(file)
            viewerClassNames.append(className)


    errors = []

    for className in viewerClassNames:
        try:
            className = string.lower(className[0]) + className[1:]
            viewerModule = imp.load_module(className, *imp.find_module(className, __path__))
            
            return viewerModule.make(vars = vars, title = title, limits = limits)
        except Exception, s:
            errors.append("%s: %s" % (className, s))
            
    raise ImportError, "Failed to import a viewer: %s" % str(errors)
