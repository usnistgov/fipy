#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "pyxviewer.py"
 #                                    created: 6/25/04 {3:17:21 PM} 
 #                                last update: 6/25/04 {3:17:21 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #  Author: Alexander Mont
 #    mail: NIST
 #     www: http://ctcms.nist.gov/
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
 # Update# Author Reason
 # -------  ---   ---------
 # 1.0      ADM   Original
 #
 # ###################################################################
 ##

"""
The PyxViewer module takes a variable as input and plots a contour map of the values of that variable. The contour map consists
of a large number of small rectangles (hereafter referred to as 'plot cells') each of which is a different color corresponding
to the value of the variable at the center of it. The PyxViewer has a subclass Grid2DPyxViewer. Grid2DPyxViewer should be used
when the variable's mesh is a Grid2D mesh, and PyxViewer should be used otherwise.

NOTE: The current implementation of PyxViewer (but not Grid2DPyxViewer) is known to cause a segmentation fault when the mesh
size and/or the number of plot cells is large.

A PyxViewer is created with one argument: the variable to be plotted. To actually plot the variable, it is necessary to
call the PyxViewer's plot() method. The plot() method has the following keyword arguments:

debug - When set to 1, prints debugging information. Default value is 0.
returnlist - When set to 1, returns the displaylist. Used for debugging and test cases. Default value is 0.
minx, maxx - The minimum and maximum values of X to appear on the plot. Default values are 0.0 and 10.0.
miny, maxy - The minimum and maximum values of Y to appear on the plot. Default values are 0.0 and 10.0.
minval, maxval - The minimum and maximum values on the color scale. A value of minval will appear completely
blue and a value of maxval will appear completely red. If no values are specified, minval will default to the lowest value of the variable
over the area to be plotted and maxval will default to the highest value of the variable over the area to be plotted/ 
resolution - The distance between the centers of adjacent plot cells. Lower values will result in higher quality images. If no value
is specified, the resolution will be set such that there will be approximately 1,000 plot cells.
filename - The name of the file to store the resulting image in. An extension is not needed here - the .eps extension
will automatically be added.
viewcmd - The OS command used to access your graphics viewer. If specified, the graphics viewer will be automatically opened
to view your plot. If no viewcmd is specified, no graphics program will open.
xlabel - The label to use for the X axis. Default is 'X values'.
ylabel - The label to use for the Y axis. Default is 'Y values'.
gridcorrect - The amount by which adjacent plot cells should overlap. Used for correcting the 'grid effect'. Default value is 0.03.

Test cases:

   >>> from fipy.meshes.grid2D import Grid2D
   >>> import fipy.variables.cellVariable
   >>> mesh = Grid2D(1.0, 1.0, 2, 2)
   >>> myvar = fipy.variables.cellVariable.CellVariable(mesh)
   >>> myvar[0] = 1.0
   >>> myvar[1] = 2.0
   >>> myvar[2] = 3.0
   >>> myvar[3] = 4.0
   >>> myviewer = Grid2DPyxViewer(myvar)
   >>> array = myviewer.getValueMatrix(0.0, 2.0, 0.0, 2.0, 1.0)
   >>> testlist = array.tolist() 
   >>> print testlist
   [[1.0, 3.0], [2.0, 4.0]]
   
   >>> from fipy.meshes.grid2D import Grid2D
   >>> import fipy.variables.cellVariable
   >>> mesh = Grid2D(1.0, 1.0, 2, 2)
   >>> myvar = fipy.variables.cellVariable.CellVariable(mesh)
   >>> myvar[0] = 1.0
   >>> myvar[1] = 2.0
   >>> myvar[2] = 3.0
   >>> myvar[3] = 4.0
   >>> myviewer = Grid2DPyxViewer(myvar)
   >>> testlist = myviewer.plot(minx = 0.0, maxx = 2.0, miny = 0.0, maxy = 2.0, resolution = 1.0, minval = 1.0, maxval = 5.0, gridcorrect = 0.0, returnlist = 1) 
   >>> print testlist
   [[0.0, 1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 1.0, 2.0, 0.5], [1.0, 2.0, 0.0, 1.0, 0.25], [1.0, 2.0, 1.0, 2.0, 0.75]]

   >>> from fipy.meshes.grid2D import Grid2D
   >>> import fipy.variables.cellVariable
   >>> mesh = Grid2D(1.0, 1.0, 2, 2)
   >>> myvar = fipy.variables.cellVariable.CellVariable(mesh)
   >>> myvar[0] = 2.0
   >>> myvar[1] = 3.0
   >>> myvar[2] = 5.0
   >>> myvar[3] = 6.0
   >>> myviewer = Grid2DPyxViewer(myvar)
   >>> testlist = myviewer.plot(minx = 0.0, maxx = 2.0, miny = 0.0, maxy = 2.0, resolution = 1.0, gridcorrect = 0.0, returnlist = 1) 
   >>> print testlist
   [[0.0, 1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 1.0, 2.0, 0.75], [1.0, 2.0, 0.0, 1.0, 0.25], [1.0, 2.0, 1.0, 2.0, 1.0]]

   >>> from fipy.meshes.grid2D import Grid2D
   >>> import fipy.variables.cellVariable
   >>> mesh = Grid2D(1.0, 1.0, 2, 2)
   >>> myvar = fipy.variables.cellVariable.CellVariable(mesh)
   >>> myvar[0] = 1.0
   >>> myvar[1] = 2.0
   >>> myvar[2] = 3.0
   >>> myvar[3] = 4.0
   >>> myviewer = Grid2DPyxViewer(myvar)
   >>> array = myviewer.getValueMatrix(0.0, 2.0, 0.0, 2.0, 0.5)
   >>> testlist = array.tolist()
   >>> print testlist
   [[0.25, 1.25, 2.25, 3.25], [0.75, 1.75, 2.75, 3.75], [1.25, 2.25, 3.25, 4.25], [1.75, 2.75, 3.75, 4.75]]

   >>> from fipy.meshes.grid2D import Grid2D
   >>> import fipy.variables.cellVariable
   >>> mesh = Grid2D(1.0, 1.0, 2, 2)
   >>> myvar = fipy.variables.cellVariable.CellVariable(mesh)
   >>> myvar[0] = 1.0
   >>> myvar[1] = 2.0
   >>> myvar[2] = 3.0
   >>> myvar[3] = 4.0
   >>> myviewer = PyxViewer(myvar)
   >>> array = myviewer.getValueMatrix(0.0, 2.0, 0.0, 2.0, 0.5)
   >>> testlist = array.tolist()
   >>> print testlist
   [[1.0, 1.0, 3.0, 3.0], [1.0, 1.0, 3.0, 3.0], [2.0, 2.0, 4.0, 4.0], [2.0, 2.0, 4.0, 4.0]]
   

   
"""

__docformat__ = 'restructuredtext'

import Numeric
import os
import pyx
import pyx.color
import time

pyx.text.set(fontmaps="psfonts.cmz")

class PyxViewer:
    
    def __init__(self, variable):
        self.var = variable
        self.theArray = Numeric.array(self.var)    ## the array of variable values, indexed by cell ID number

## ------------------------------------------------------------------------------------

    def plot(self, returnlist = 0, debug=0, minx=0.0, maxx=10.0, miny=0.0, maxy=10.0, minval=None, maxval=None, resolution = None, filename = None, viewcmd = None, xlabel = "X values", ylabel = "Y values", gridcorrect = 0.03):
        
        ## initialize variables
        
        starttime = time.time()
        thepalette = pyx.color.palette.ReverseRainbow
        ## calculate the resolution. If there is no resolution set, calculate a "default" resolution that will lead to 1,000 points being plotted.
        ## If a resolution is given, use that.
        
        if(resolution == None):
            resolution = ((((maxx - minx)*(maxy - miny))/1000.0) ** 0.5)

        ## calculate the number of points to plot in the X direction (xsize) and the Y direction (ysize). The total number of points plotted
        ## will be equal to (xsize * ysize).
        xsize = ((maxx - minx) / resolution) - 0.5
        ysize = ((maxx - minx) / resolution) - 0.5
        if (xsize == int(xsize)):
            xsize = int(xsize)
        else:
            xsize = int(xsize) + 1
        if (ysize == int(ysize)):
            ysize = int(ysize)
        else:
            ysize = int(ysize) + 1            

        ## xarr, yarr, and valarr are all arrays with the dimensions (xsize, ysize). Each element in each array represents one point to be plotted.
        ## Elements in xarr represent the X coordinates of the points, those in yarr represent the Y coordinates of the points, and those in
        ## valarr represent the values at those points.
        xarr = Numeric.fromfunction(lambda x, y: minx + ((x+0.5) * resolution), (xsize, ysize))
        yarr = Numeric.fromfunction(lambda x, y: miny + ((y+0.5) * resolution), (xsize, ysize)) 
        valarr = self.getValueMatrix(minx, maxx, miny, maxy, resolution)
        
        ## This section takes the values above and generates the list to put into the pyx.graph.plot method.
        ## If maxval and minval are not given, generate them.
        xlist = Numeric.reshape(xarr, (Numeric.size(xarr),))
        ylist = Numeric.reshape(yarr, (Numeric.size(yarr),))
        vallist = Numeric.reshape(valarr, (Numeric.size(valarr),))
        if(minval == None):
            minval = min(vallist)
        if(maxval == None):
            maxval = max(vallist)
        vallist = (vallist - minval) / (maxval - minval)
        minxlist= xlist - (0.5 * resolution)
        maxxlist = xlist + (0.5 * resolution)
        minylist = ylist - (0.5 * resolution)
        maxylist = ylist + (0.5 * resolution)
        baselist = Numeric.array([minxlist, (maxxlist + gridcorrect), minylist, (maxylist + gridcorrect), vallist])
        displaylist = Numeric.transpose(baselist)

        ## If the user did not specify a filename but does want to view the picture, create a file anyway so the picture can be viewed.
        if(filename == None):
            if(viewcmd != None):
                filename = "temp"
        displayinput = displaylist.tolist()

        ## Create the graph and plot it.
        if(filename != None):
            mygraph = pyx.graph.graphxy(height = 8, width = 8, x = pyx.graph.axis.linear(min = minx, max = maxx, title = xlabel), y = pyx.graph.axis.linear(min = miny, max = maxy, title = ylabel))
            mygraph.plot(pyx.graph.data.list(displayinput, xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5), pyx.graph.style.rect(pyx.color.palette.Rainbow))
            mygraph.dodata()
            mygraph.dokey()
            mygraph.writeEPSfile(filename)
        tottime = time.time() - starttime
        if(debug == 1):
            print("time to execute:")
            print tottime
            print("resolution:")
            print resolution

        if(viewcmd != None):
            os.system(viewcmd + " " + filename + ".eps")

        if(returnlist == 1):
            return displayinput
        else:
            return None
## ---------------------------------------------------------------------------------------------------
            
    def getValueMatrix(self, minx, maxx, miny, maxy, resolution):
        xsize = ((maxx - minx) / resolution) - 0.5
        ysize = ((maxx - minx) / resolution) - 0.5
        if (xsize == int(xsize)):
            xsize = int(xsize)
        else:
            xsize = int(xsize) + 1
        if (ysize == int(ysize)):
            ysize = int(ysize)
        else:
            ysize = int(ysize) + 1    
        valarr = Numeric.fromfunction(lambda x, y: self.getValue(minx + ((x+0.5) * resolution), miny + ((y+0.5) * resolution)), (xsize, ysize))
        return valarr
    
## ---------------------------------------------------------------------------------------------------
            
    def getValue(self, x, y):

        ## This function gets the value of the variable at a given point x, y. This function is designed to be used with
        ## x and y being arrays.
        
        cellIDarray = self.var.getMesh().getNearestCellID((x, y))
        resArray = Numeric.take(self.theArray, cellIDarray)
        return resArray

##------------------##------------------##------------------------------## ----------------##----------------##----------

InvalidClassException = "Wrong class"

class Grid2DPyxViewer(PyxViewer):
    def __init__(self, variable):
        if(variable.getMesh().__class__.__name__ != "Grid2D"):
            raise InvalidClassException, "Grid2DPyxViewer can be used only with Grid2D meshes"
        else:
            self.var = variable
            self.theArray = Numeric.array(self.var)

    def getValue(self, x, y):
        ## note - need to change so it takes into account that the lower left corner might not be (0,0)
        ## This function uses the method of bilinear interpolation on the grid square. This method is described in
        ## "Numerical Recipes in C", William H. Press et al., page 105.
        width = self.var.getMesh().__getstate__()["nx"]
        height = self.var.getMesh().__getstate__()["ny"] 
        lowerleft = self.var.getMesh().getCells()[0].getCenter()
        startx = lowerleft[0]
        starty = lowerleft[1]
        dx = self.var.getMesh().getMeshSpacing()[0]
        dy = self.var.getMesh().getMeshSpacing()[1]
        xcells = ((x - startx) / dx)
        ycells = ((y - starty) / dy)
        onLeftSide = (xcells < 0)
        onRightSide = (xcells >= (width - 1.0))
        onBottom = (ycells < 0)
        onTop = (ycells >= (height - 1.0))
        xnumcells = xcells.astype(Numeric.Int)
        ynumcells = ycells.astype(Numeric.Int)
        xnumcells = xnumcells - onRightSide
        ynumcells = ynumcells - onTop
        Tvalues = (xcells - xnumcells)
        Uvalues = (ycells - ynumcells)
        Y1coeffs = (1 - Tvalues)*(1 - Uvalues)
        Y2coeffs = Tvalues * (1 - Uvalues)
        Y3coeffs = Tvalues * Uvalues
        Y4coeffs = (1 - Tvalues) * Uvalues
        cellIDarray = (ynumcells * width) + xnumcells
        resArray = Numeric.take(self.theArray, cellIDarray) * Y1coeffs   ## contribution from lower left point
        resArray = resArray + Numeric.take(self.theArray, cellIDarray + 1) * Y2coeffs   ## contribution from lower right point
        resArray = resArray + Numeric.take(self.theArray, cellIDarray + (width + 1)) * Y3coeffs   ## contribution from upper right point
        resArray = resArray + Numeric.take(self.theArray, cellIDarray + width) * Y4coeffs   ##contribution from upper left point
        return resArray


def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
