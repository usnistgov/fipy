#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gist2DViewer.py"
 #                                    created: 11/10/03 {2:48:25 PM} 
 #                                last update: 2/21/07 {12:24:17 PM} 
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.viewers.gistViewer.gistViewer import GistViewer

class Gist2DViewer(GistViewer):
    """
    Displays a contour plot of a 2D `CellVariable` object.

    """
    
    def __init__(self, vars, limits = None, title = None, palette = 'heat.gp', grid = 1, dpi = 75):
        """
        Creates a `Gist2DViewer`.
        
        >>> from fipy import *
        >>> mesh = Grid2D(nx=50, ny=100, dx=0.1, dy=0.01)
        >>> x, y = mesh.getCellCenters()[...,0], mesh.getCellCenters()[...,1]
        >>> var = CellVariable(mesh=mesh, name=r"$sin(x y)$", value=numerix.sin(x * y))
        >>> vw = Gist2DViewer(vars=var, 
        ...                   limits={'ymin':10, 'ymax':90, 'datamin':-0.9, 'datamax':2.0},
        ...                   title="Gist2DViewer test")
        >>> if locals().has_key('vw'):
        ...     vw.plot()
        ...     raw_input("Describe any problems with this figure or hit Return: ").strip()
        ...     del vw
        ''

        >>> mesh = Tri2D(nx=50, ny=100, dx=0.1, dy=0.01)
        >>> x, y = mesh.getCellCenters()[...,0], mesh.getCellCenters()[...,1]
        >>> var = CellVariable(mesh=mesh, name=r"$sin(x y)$", value=numerix.sin(x * y))
        >>> vw = Gist2DViewer(vars=var, 
        ...                   limits={'ymin':10, 'ymax':90, 'datamin':-0.9, 'datamax':2.0},
        ...                   title="Gist2DViewer test")
        >>> if locals().has_key('vw'):
        ...     vw.plot()
        ...     raw_input("Describe any problems with this figure or hit Return: ").strip()
        ...     del vw

        :Parameters:
          - `vars`: A `CellVariable` or tuple of `CellVariable` objects to plot.
            Only the first 2D `CellVariable` will be plotted.
          - `limits`: A dictionary with possible keys `'xmin'`, `'xmax'`, 
            `'ymin'`, `'ymax'`, `'datamin'`, `'datamax'`. Any limit set to 
            a (default) value of `None` will autoscale.
          - `title`: Displayed at the top of the Viewer window.
          - `palette`: The color scheme to use for the image plot. Default is 
            `heat.gp`. Another choice would be `rainbow.gp`.
          - `grid`: Whether to show the grid lines in the plot. Default is 1. 
            Use 0 to switch them off.
            
        """
        GistViewer.__init__(self, vars = vars, limits = limits, 
                            title = " ", dpi = dpi)
                            
        self.palette = palette
        self.grid = grid
        
    def _getSuitableVars(self, vars):
        from fipy.variables.cellVariable import CellVariable
        vars = [var for var in GistViewer._getSuitableVars(self, vars) \
          if (var.getMesh().getDim() == 2 and isinstance(var, CellVariable))]
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
            gist.limits(self._getLimit('xmin'), self._getLimit('xmax'), self._getLimit('ymin'), self._getLimit('ymax'))

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

        
        vertexIDs = self.mesh._getOrderedCellVertexIDs().flat

        import MA
        
        if type(vertexIDs) is type(MA.array(0)):
            vertexIDs = vertexIDs.compressed()
            
        vertexCoords = self.mesh.getVertexCoords()

        xCoords = numerix.take(vertexCoords[:,0], numerix.array(vertexIDs).flat)
        yCoords = numerix.take(vertexCoords[:,1], numerix.array(vertexIDs).flat)

        import gist

        import Numeric
        gist.plfp(Numeric.array(numerix.array(self.vars[0])), yCoords, xCoords, self.mesh._getNumberOfFacesPerCell(), cmin=datamin, cmax=datamax)

        import colorbar

        colorbar._color_bar(minz=datamin, maxz=datamax, ncol=240, zlabel=self.vars[0].getName())

        GistViewer.plot(self, filename = filename)

    def plotMesh(self, filename = None):
        """
        Plot the `CellVariable`'s mesh as a wire frame.
        """
        self._plot()
        
        faceVertexIDs = self.mesh._getFaceVertexIDs()
        vertexCoords = self.mesh.getVertexCoords()
        
        x0 = numerix.take(vertexCoords[:,0], faceVertexIDs[:,0])
        y0 = numerix.take(vertexCoords[:,1], faceVertexIDs[:,0])
        x1 = numerix.take(vertexCoords[:,0], faceVertexIDs[:,1])
        y1 = numerix.take(vertexCoords[:,1], faceVertexIDs[:,1])
        
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

