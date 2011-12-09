#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gist1DViewer.py"
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

from fipy.tools import numerix
from fipy.viewers.gistViewer.gistViewer import _GistViewer
from fipy.tools.decorators import deprecateGist
__all__ = ["Gist1DViewer"]

@deprecateGist
class Gist1DViewer(_GistViewer):
    """Displays a y vs. x plot of one or more 1D `CellVariable` objects.
    """
    
    __doc__ += _GistViewer._test1D(viewer="Gist1DViewer")

    def __init__(self, vars, title=None, xlog=0, ylog=0, style="work.gs", limits={}, **kwlimits):
        """
        Creates a `Gist1DViewer`.
        
        :Parameters:
          vars
            a `CellVariable` or tuple of `CellVariable` objects to plot
          title
            displayed at the top of the `Viewer` window
          xlog
            log scaling of x axis if `True`
          ylog
            log scaling of y axis if `True`
          stye
            the Gist stylefile to use.
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          xmin, xmax, datamin, datamax
            displayed range of data. Any limit set to 
            a (default) value of `None` will autoscale.
            (*ymin* and *ymax* are synonyms for *datamin* and *datamax*).
        """
        kwlimits.update(limits)
        _GistViewer.__init__(self, vars=vars, title=title, **kwlimits)
        
        self.xlog = xlog
        self.ylog = ylog
        self.style = style
        
    def _getSuitableVars(self, vars):
        from fipy.variables.cellVariable import CellVariable
        vars = [var for var in _GistViewer._getSuitableVars(self, vars) \
          if (var.mesh.dim == 1 and isinstance(var, CellVariable))]
        if len(vars) > 1:
            vars = [var for var in vars if var.mesh is vars[0].mesh]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError, "Can only plot 1D data"
        return vars
        
    def _getArrays(self):
        arrays = []
        
        for var in self.vars:
            arrays.append((numerix.array(var), numerix.array(var.mesh.cellCenters[0])))
            
        return arrays
        
    def _plotArrays(self):
        import gist
        
        for array in self._getArrays():
            gist.plg(*array)
            
        gist.logxy(self.xlog, self.ylog)

    def plot(self, filename = None):
        import gist

        gist.window(self.id, wait = 1, style = self.style)
        gist.pltitle(self.title)
        gist.animate(1)

        if self.limits != None:
            gist.limits(self._getLimit('xmin'), 
                        self._getLimit('xmax'), 
                        self._getLimit(('datamin', 'ymin')), 
                        self._getLimit(('datamax', 'ymax')))
            
        self._plotArrays()
            
        _GistViewer.plot(self, filename = filename)

if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
