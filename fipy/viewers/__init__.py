import os
import sys
import string
import glob
import imp

# what about vector variables?

def make(vars, title = None, limits = None):
    if type(vars) not in [type([]), type(())]:
        vars = [vars]
    
##     if os.environ.has_key('FiPy_%dD_ViewerClass' % dim):
##         viewerClassNames = [os.environ['FiPy_%dD_ViewerClass' % dim]]
##     else:
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
