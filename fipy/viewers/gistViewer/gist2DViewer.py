#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gist2DViewer.py"
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
__all__ = ["Gist2DViewer"]

@deprecateGist
class Gist2DViewer(_GistViewer):
    """Displays a contour plot of a 2D `CellVariable` object.
    """
    
    __doc__ += _GistViewer._test2D(viewer="Gist2DViewer")
    __doc__ += _GistViewer._test2Dirregular(viewer="Gist2DViewer")
    
    def __init__(self, vars, title=None, palette='heat.gp', grid=True, dpi=75, limits={}, **kwlimits):
        """Creates a `Gist2DViewer`.
        
        :Parameters:
          vars
            a `CellVariable` object.
          title
            displayed at the top of the `Viewer` window
          palette
            the color scheme to use for the image plot. Default is 
            `heat.gp`. Another choice would be `rainbow.gp`.
          grid
            whether to show the grid lines in the plot.
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          xmin, xmax, ymin, ymax, datamin, datamax
            displayed range of data. Any limit set to 
            a (default) value of `None` will autoscale.
        """
        kwlimits.update(limits)
        _GistViewer.__init__(self, vars=vars, title=" ", dpi=dpi, **kwlimits)
                            
        self.palette = palette
        self.grid = grid
        
    def _getSuitableVars(self, vars):
        from fipy.variables.cellVariable import CellVariable
        vars = [var for var in _GistViewer._getSuitableVars(self, vars) \
          if (var.mesh.dim == 2 and isinstance(var, CellVariable))]
        if len(vars) == 0:
            from fipy.viewers import MeshDimensionError
            raise MeshDimensionError, "Can only plot 2D data"
        # this viewer can only display one variable
        return [vars[0]]

    def _plot(self):
        import gist
    
        gist.window(self.id, wait = 1)
        gist.animate(1)
        gist.pltitle(self.title)
        gist.palette(self.palette)
        gist.gridxy(self.grid)

        if self.limits != None:
            gist.limits(self._getLimit('xmin'), 
                        self._getLimit('xmax'), 
                        self._getLimit('ymin'), 
                        self._getLimit('ymax'))

    def plot(self, filename = None):
        """
        Plot the `CellVariable` as a contour plot.
        """
        self._plot()

        datamin = self._getLimit(('datamin', 'zmin'))
        datamax = self._getLimit(('datamax', 'zmax'))
        
        if datamin == 'e':
            datamin = None
        
        if datamax == 'e':
            datamax = None
            
        datamin, datamax = self._autoscale(vars=self.vars,
                                           datamin=datamin,
                                           datamax=datamax)

        if datamax == datamin:
            datamax = datamin + 1e-10

        vertexIDs = self.mesh._orderedCellVertexIDs

        vertexCoords = self.mesh.vertexCoords

        xCoords = numerix.take(vertexCoords[0], vertexIDs).flatten("FORTRAN")
        yCoords = numerix.take(vertexCoords[1], vertexIDs).flatten("FORTRAN")

        import gist

        import Numeric
        gist.plfp(Numeric.array(numerix.array(self.vars[0])), yCoords, xCoords, self.mesh._numberOfFacesPerCell, cmin=datamin, cmax=datamax)

        import fipy.viewers.gistViewer.colorbar

        colorbar._color_bar(minz=datamin, maxz=datamax, ncol=240, zlabel=self.vars[0].name)

        _GistViewer.plot(self, filename = filename)

    def plotMesh(self, filename = None):
        self._plot()
        
        faceVertexIDs = self.mesh.faceVertexIDs
        vertexCoords = self.mesh.vertexCoords
        
        x0 = numerix.take(vertexCoords[0], faceVertexIDs[0])
        y0 = numerix.take(vertexCoords[1], faceVertexIDs[0])
        x1 = numerix.take(vertexCoords[0], faceVertexIDs[1])
        y1 = numerix.take(vertexCoords[1], faceVertexIDs[1])
        
        import gist
        
        gist.pldj(x0, y0, x1, y1)

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

if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
