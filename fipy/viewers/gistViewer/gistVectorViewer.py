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

from fipy.viewers.gistViewer.gistViewer import _GistViewer
from fipy.variables.cellVariable import CellVariable
from fipy.variables.faceVariable import FaceVariable
from fipy.tools import numerix
from fipy.tools.decorators import deprecateGist
__all__ = ["GistVectorViewer"]

@deprecateGist
class GistVectorViewer(_GistViewer):
    """Displays a vector plot of a 2D rank-1 `CellVariable` or
    `FaceVariable` object using gist.
    """

    __doc__ += _GistViewer._test2Dvector(viewer="GistVectorViewer")
    __doc__ += _GistViewer._test2DvectorIrregular(viewer="GistVectorViewer")
    
    def __init__(self, vars, title=None, limits={}, **kwlimits):
        """Creates a `GistVectorViewer`.
        
        :Parameters:
          vars
            a rank-1 `CellVariable` or `FaceVariable` object.
          title
            displayed at the top of the `Viewer` window
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. Any limit set to 
            a (default) value of `None` will autoscale.
        """
        kwlimits.update(limits)
	_GistViewer.__init__(self, vars=vars, title=title, **kwlimits)
        
    def _getSuitableVars(self, vars):
        vars = [var for var in _GistViewer._getSuitableVars(self, vars) \
          if (var.mesh.dim == 2 \
              and (isinstance(var, FaceVariable) \
                   or isinstance(var, CellVariable)) and var.rank == 1)]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError, "Can only plot 2D vector data"
        # this viewer can only display one variable
        return [vars[0]]
        
    def plot(self, filename=None):
        import gist

        gist.window(self.id, wait = 1)
        gist.pltitle(self.title)
        gist.animate(1)
        
        var = self.vars[0]
        
        if isinstance(var, FaceVariable):
            x, y = var.mesh.faceCenters
        elif isinstance(var, CellVariable):
            x, y = var.mesh.cellCenters
        
        gist.plmesh(numerix.array([y, y]), numerix.array([x, y]))

        vx = numerix.array(var[0])
        vy = numerix.array(var[1])
        
        maxVec = var.mag.max().value
        maxGrid = var.mesh._cellDistances.max()
        
        gist.plv(numerix.array([vy,vy]), numerix.array([vx,vx]), scale=maxGrid / maxVec * 3, hollow=1, aspect=0.25) #,scale=0.002)
        
        if filename is not None:
            
            gist.hcp_file(filename)
            gist.hcp()

        gist.fma()
        
    def getArray(self):
        pass
        
if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
