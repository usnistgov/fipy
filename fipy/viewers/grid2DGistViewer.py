#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gistViewer.py"
 #                                    created: 11/10/03 {2:48:25 PM} 
 #                                last update: 9/1/04 {6:33:57 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #  Author: James Warren
 #  E-mail: jwarren@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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

The `Grid2DGistViewer` uses the gist plotting package to display a
contour plot of a variable. It assumes the mesh that the variable has
is a `Grid2D` type mesh.

It takes a resolution argument when it is instantiated. This argument increases the
resolution of the plot by doing interpolation on the values from the variable.
Here is a test for increasing the resolution:

   >>> from fipy.variables.cellVariable import CellVariable
   >>> from fipy.meshes.grid2D import Grid2D
   >>> nx = 2
   >>> ny = 4
   >>> dx = .5
   >>> dy = 1.
   >>> mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)
   >>> v = (mesh.getCellCenters()[:,0] + mesh.getCellCenters()[:,1])
   >>> var = CellVariable(mesh = mesh, value = v)
   >>> resolution = 3
   >>> viewer = Grid2DGistViewer(var = var, resolution = resolution)
   >>> pow = 2**(resolution-1)
   >>> nx = pow * nx - (pow - 1)
   >>> ny = pow * ny - (pow - 1)
   >>> index = Numeric.indices((ny, nx))
   >>> ndx = dx / pow
   >>> ndy = dy / pow
   >>> answer = (index[1] * ndx + dx / 2) + (index[0] * ndy + dy / 2)
   >>> Numeric.allclose(answer, viewer.getArray())
   1
   
"""

import Numeric
 
from fipy.viewers.gistViewer import GistViewer

class Grid2DGistViewer(GistViewer):
    
    def __init__(self, var = None, minVal=None, maxVal=None, palette = 'heat.gp', grid = 1, resolution = 1, limits = None, dpi = 75):
        """
        The following arguments can be given to a `Grid2DGistViewer`:

        `var` - The variable that is to be plotted

        `minVal` - The minimum value to display in the plot. The
        default will use the minimum value from the variable values.

        `maxVal` - The maximum value to display in the plot. The
        default will use the maximum value from the variable values.

        `palette` - The color scheme to use for the contour
        plot. Default is `heat.gp`. Another choice would be
        `rainbow.gp`.

        `grid` - Whether to show the grid lines in the plot. Default
        is 1. Use 0 to switch them off.

        `resolution` - Is an integer greater than 0. It corresponds to
        the amount of resolution required for the plot. The default is
        1. The increase in resolution is equivalent to `2**(resolution - 1)`.

        `limits` - Plot limits of interest, tuple of four values
        corresponding to `(xmin, xmax, ymin, ymax)`
        
        """

        self.var = var

        self.resolution = resolution
        GistViewer.__init__(self, minVal, maxVal, title = self.var.name, palette = palette, grid = grid, limits = limits, dpi = dpi)

    def setVar(self,var):
        self.var = var

    def getArray(self):

        nx,ny = self.var.getMesh().getShape()
        array = Numeric.reshape(self.var.getNumericValue(),(ny,nx))

        for i in range(self.resolution - 1):
            array = self.increaseResolution(array)

        return array

    def increaseResolution(self, array):
        """
        >>> from fipy.variables.cellVariable import CellVariable
        >>> from fipy.meshes.grid2D import Grid2D

        >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 4)
        >>> v = (mesh.getCellCenters()[:,0] + mesh.getCellCenters()[:,1])
        >>> var = CellVariable(mesh = mesh, value = v)
        >>> viewer = Grid2DGistViewer(var = var)
        >>> nx,ny = mesh.getShape()
        >>> array = Numeric.reshape(var.getNumericValue(),(ny,nx))
        >>> print viewer.increaseResolution(array)
        [[ 1. , 1.5, 2. ,]
         [ 1.5, 2. , 2.5,]
         [ 2. , 2.5, 3. ,]
         [ 2.5, 3. , 3.5,]
         [ 3. , 3.5, 4. ,]
         [ 3.5, 4. , 4.5,]
         [ 4. , 4.5, 5. ,]]

        """

        return self._increaseResolution(self._increaseResolution(array), axis = 1)

    def _increaseResolution(self, array, axis = 0):
        """
        >>> from fipy.variables.cellVariable import CellVariable
        >>> from fipy.meshes.grid2D import Grid2D

        >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 4)
        >>> v = (mesh.getCellCenters()[:,0] + mesh.getCellCenters()[:,1])
        >>> var = CellVariable(mesh = mesh, value = v)
        >>> viewer = Grid2DGistViewer(var = var)
        >>> nx,ny = mesh.getShape()
        >>> array = Numeric.reshape(var.getNumericValue(),(ny,nx))
        >>> print viewer._increaseResolution(array)
        [[ 1. , 2. ,]
         [ 1.5, 2.5,]
         [ 2. , 3. ,]
         [ 2.5, 3.5,]
         [ 3. , 4. ,]
         [ 3.5, 4.5,]
         [ 4. , 5. ,]]
        >>> print viewer._increaseResolution(array, axis = 1)
        [[ 1. , 1.5, 2. ,]
         [ 2. , 2.5, 3. ,]
         [ 3. , 3.5, 4. ,]
         [ 4. , 4.5, 5. ,]]

        """
        
        n = array.shape[axis]

        if axis == 0:
            interpolatedArray = (array[1:,:] + array[:-1,:]) / 2
        else:
            interpolatedArray = (array[:,1:] + array[:,:-1]) / 2

        array = Numeric.concatenate((array, interpolatedArray), axis = axis)

        ind = Numeric.zeros((n, 2))
        ind[:,0] = Numeric.arange(n)
        ind[:,1] = Numeric.arange(n) + n
        ind = ind.flat[:-1]

        return Numeric.take(array, ind, axis = axis)

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
