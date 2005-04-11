#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gist2DViewer.py"
 #                                    created: 11/10/03 {2:48:25 PM} 
 #                                last update: 4/5/05 {5:33:37 PM} 
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

import Numeric
import fipy.tools.array 
from fipy.viewers.gistViewer import GistViewer

import gist
import colorbar

class Gist2DViewer(GistViewer):
    """
    Displays a contour plot of a 2D `CellVariable` object.
    Usage:

    ::
    
        viewer = Gist2DViewer(var)
        viewer.plot()

    """
    
    def __init__(self, vars, limits = None, title = None, palette = 'heat.gp', grid = 1, dpi = 75):
        """
        Creates a `Gist2DViewer`.
        
        :Parameters:
          - `vars` : A `CellVariable` object to plot (not a list or tuple).
          - `limits` : A dictionary with possible keys `xmin`, `xmax`, 
            `ymin`, `ymax`, `zmin`, `zmax`, `datamin`, `datamax`.  A 1D
            Viewer will only use `xmin` and `xmax`, a 2D viewer will also
            use `ymin` and `ymax`, and so on.  All viewers will use
            `datamin` and `datamax`.  Any limit set to a (default) value of
            `None` will autoscale.
          - `title` : Displayed at the top of the Viewer window.
          - `palette` : The color scheme to use for the image plot. Default is 
            `heat.gp`. Another choice would be `rainbow.gp`.
          - `grid` : Whether to show the grid lines in the plot. Default is 1. 
            Use 0 to switch them off.
            
        """

        GistViewer.__init__(self, vars = vars, limits = limits, title = title, dpi = dpi)

        if len(self.vars) != 1:
            raise IndexError, "A 2D Gist viewer can only display one CellVariable"
        
        self.palette = palette
        self.grid = grid

    def _plot(self):
        gist.window(self.id, wait = 1)
        gist.animate(1)
        gist.pltitle(self.title)
        gist.palette(self.palette)
        gist.gridxy(self.grid)

        if self.limits != None:
            gist.limits(self._getLimit('xmin'), self._getLimit('xmax'), self._getLimit('ymin'), self._getLimit('ymax'))

    def plot(self):
        """
        Plot the `CellVariable` as a contour plot.
        """
        self._plot()

        minVal = self._getLimit('datamin')
        maxVal = self._getLimit('datamax')
        
        if minVal == 'e':
            minVal = fipy.tools.array.min(self.vars[0][:])
            for var in self.vars[1:]:
                minVal = min(minVal, fipy.tools.array.min(var[:]))

        if maxVal == 'e':
            maxVal = fipy.tools.array.max(self.vars[0][:])
            for var in self.vars[1:]:
                maxVal = max(maxVal, fipy.tools.array.max(var[:]))

        if maxVal == minVal:
            maxVal = minVal + 1e-10
            
        vertexIDs = self.mesh._getOrderedCellVertexIDs()
        Nfac = self.mesh._getMaxFacesPerCell()
        Ncells = self.mesh.getNumberOfCells()
        xCoords = Numeric.take(self.mesh.getVertexCoords()[:,0], vertexIDs.flat)
        yCoords = Numeric.take(self.mesh.getVertexCoords()[:,1], vertexIDs.flat)
        gist.plfp(Numeric.array(self.vars[0]), yCoords, xCoords, Nfac * Numeric.ones(Ncells), cmin = minVal, cmax = maxVal)

        colorbar._color_bar(minz = minVal, maxz = maxVal, ncol=240, zlabel = 'fred')

        gist.fma()

    def plotMesh(self):
        """
        Plot the `CellVariable`'s mesh as a wire frame.
        """
        self._plot()
        
        faceVertexIDs = self.mesh._getFaceVertexIDs()
        vertexCoords = self.mesh.getVertexCoords()
        
        from fipy.tools import array
        
        x0 = array.take(vertexCoords[:,0], faceVertexIDs[:,0])
        y0 = array.take(vertexCoords[:,1], faceVertexIDs[:,0])
        x1 = array.take(vertexCoords[:,0], faceVertexIDs[:,1])
        y1 = array.take(vertexCoords[:,1], faceVertexIDs[:,1])
        
        gist.pldj(x0, y0, x1, y1)

        gist.fma()

