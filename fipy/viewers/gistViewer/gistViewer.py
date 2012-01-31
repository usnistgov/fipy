#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gistViewer.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.viewers.viewer import AbstractViewer

class _GistViewer(AbstractViewer):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    
    _id=0
    
    def __init__(self, vars, title=None, dpi=75, **kwlimits):
        """
        Create a `_GistViewer` object.
        
        :Parameters:
          vars
            a `CellVariable` or tuple of `CellVariable` objects to plot
          title
            displayed at the top of the `Viewer` window
          dpi
            the dot-per-inch resolution of the display
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. A 1D `Viewer` will only use `xmin` and
            `xmax`, a 2D viewer will also use `ymin` and `ymax`. All
            viewers will use `datamin` and `datamax`. Any limit set to a
            (default) value of `None` will autoscale.
        """
        if self.__class__ is _GistViewer:
            raise NotImplementedError, "can't instantiate abstract base class"
    
        AbstractViewer.__init__(self, vars=vars, title=title, **kwlimits)
        
        self.mesh = self.vars[0].mesh

        self.id = _GistViewer._id 
        _GistViewer._id += 1
        
        import gist
        
        gist.window(self.id, wait = 1, dpi = dpi, display = '')

    def _getLimit(self, keys, default='e'):
        return AbstractViewer._getLimit(self, keys=keys, default=default)
        
    def plot(self, filename=None):
        import gist
    
        if filename is not None:
            import os.path
            root, ext = os.path.splitext(filename)
            if ext.lower() in (".eps", ".epsi"):
                gist.eps(root)
            else:
                gist.hcp_file(filename, dump = 1)
                gist.hcp()
                gist.hcp_finish(-1)
        gist.fma()
        
    def _validFileExtensions(self):
        return [".cgm", ".eps", ".epsi"]
