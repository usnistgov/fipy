__version__ = "1.0a3"

import sys
import imp
import os
import glob

fipyPath = __path__[0]

class FiPyImporter:
    """
    We want to allow users to write
    
        >>> from fipy import Grid2D
        
    instead of
    
        >>> from fipy.meshes.grid2D import Grid2D
        
    although the latter will always work
    """

    modules = {
    }
    
    def __init__(self, path=fipyPath):
        if path != fipyPath:
            # if out class is on sys.path_hooks, we must raise
            # ImportError for any path item that we can't handle.
            raise ImportError
        self.path = path
        
    def find_module(self, fullname, path=None):
        parts = fullname.split(".")
        
        # This importer is only appropriate for short-cut finding classes of fipy.
        # We reject any of:
        #   - top-level module (path == None)
        #   - module that is not fipy (self.path not in path)
        #   - module that is buried in fipy (len(parts) != 2)
        #   - module instead of class (first character is lower-case)
        if not path or self.path not in path or len(parts) != 2 or parts[-1][0].islower():
            return None
            
        moduleName = "%s%s" % (parts[-1][0].lower(), parts[-1][1:])
        moduleFiles = ["%s%s" % (moduleName, suffix) \
            for suffix, mode, moduleType in imp.get_suffixes() \
            if moduleType is not imp.PY_COMPILED]        
        
##             if moduleType is imp.PY_SOURCE:
##                 pySuffix = suffix
##                 break
                
##         pyFile = "%s%s" % (moduleName, pySuffix)
        
        for aPath in path:
            for root, dirs, files in os.walk(aPath):
                for moduleFile in moduleFiles:
                    if moduleFile in files:
                        self.modules[fullname] = ".".join(list(os.path.split(root)) + [moduleName])
                        return self
                    
        return None

    def load_module(self, fullname):
        className = fullname.split(".")[-1]
        
        mod = __import__(self.modules[fullname], globals(), locals(), [className])
        mod = getattr(mod, className)
        
        sys.modules[fullname] = mod
        mod.__loader__ = self

        return mod
                    
sys.meta_path.append(FiPyImporter(__path__[0]))

                    
