#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh2DViewer.py"
 #                                    created: 11/10/03 {2:48:25 PM} 
 #                                last update: 10/19/04 {2:52:54 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James Warren <jwarren@nist.gov>
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

"""

This viewer displays a general unstructured 2D mesh using gist. Currently the elements must
have the same number of faces but this can be changed.
   
"""
__docformat__ = 'restructuredtext'

import Numeric
 
from fipy.viewers.gistViewer import GistViewer

try:
    import gist
except:
    print "Unable to load gist"

import fipy.tools.array

class Mesh2DGistViewer(GistViewer):
    
    def __init__(self, var = None, minVal = None, maxVal = None, palette = 'heat.gp', grid = 1, resolution = 1, limits = None, dpi = 75):
        """
	:Parameters:
	  - `var`: The variable that is to be plotted
	  - `minVal`: The minimum value to display in the plot. The
	    default will use the minimum value from the variable values.
	  - `maxVal`: The maximum value to display in the plot. The
	    default will use the maximum value from the variable values.
	  - `palette`: The color scheme to use for the contour
	    plot. Default is `heat.gp`. Another choice would be `rainbow.gp`.
	  - `grid`: Whether to show the grid lines in the plot. Default
	    is 1. Use 0 to switch them off.
	  - `limits`: Plot limits of interest. A tuple of four values corresponding 
	    to `(xmin, xmax, ymin, ymax)`
        """
        self.var = var
        GistViewer.__init__(self, minVal = minVal, maxVal = maxVal, title = self.var.name, palette = palette, grid = grid, limits = limits, dpi = dpi)

    def getArray(self):
        return Numeric.array(self.var)

    def _plot(self, array, minVal, maxVal):
        mesh = self.var.getMesh()
        vertexIDs = mesh.getOrderedCellVertexIDs()
        Nfac = mesh.getMaxFacesPerCell()
        Ncells = mesh.getNumberOfCells()
        xCoords = Numeric.take(mesh.getVertexCoords()[:,0], vertexIDs.flat)
        yCoords = Numeric.take(mesh.getVertexCoords()[:,1], vertexIDs.flat)
        gist.plfp(array, yCoords, xCoords, Nfac * Numeric.ones(Ncells), cmin = minVal, cmax = maxVal)

class Mesh2DMeshViewer(GistViewer):
    def __init__(self, mesh, minVal = None, maxVal = None, palette = 'heat.gp', grid = 1, resolution = 1, limits = None, dpi = 75):
        self.mesh = mesh
        GistViewer.__init__(self, minVal = minVal, maxVal = maxVal, title = '', palette = palette, grid = grid, limits = limits, dpi = dpi)

    def getArray(self):
        return Numeric.zeros(self.mesh.getNumberOfCells())
    
    def _plot(self, array, minVal, maxVal):
        faceVertexIDs = self.mesh.getFaceVertexIDs()
        vertexCoords = self.mesh.getVertexCoords()
        x0 = fipy.tools.array.take(vertexCoords[:,0], faceVertexIDs[:,0])
        y0 = fipy.tools.array.take(vertexCoords[:,1], faceVertexIDs[:,0])
        x1 = fipy.tools.array.take(vertexCoords[:,0], faceVertexIDs[:,1])
        y1 = fipy.tools.array.take(vertexCoords[:,1], faceVertexIDs[:,1])
        gist.pldj(x0, y0, x1, y1)
