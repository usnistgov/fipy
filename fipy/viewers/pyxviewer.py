#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "pyxviewer.py"
 #                                    created: 6/25/04 {3:17:21 PM} 
 #                                last update: 10/26/04 {10:56:28 PM} 
 #  Author: Alexander Mont <alexander.mont@nist.gov>
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

A PyxViewer is created with one argument: the variable to be plotted. To actually plot the variable, it is necessary to
call the PyxViewer's plot() method. The plot() method has the following keyword arguments:

debug - When set to 1, prints debugging information. Default value is 0.

returnlist - When set to 1, returns the displaylist. Used for debugging and test cases. Default value is 0.

minx, maxx - The minimum and maximum values of X to appear on the plot. By default, minx and maxx are set to the minimum and maximum values of the X coordinates of the mesh vertices. 

miny, maxy - The minimum and maximum values of Y to appear on the plot. By default, miny and maxy are set to the minimum and maximum values of the Y coordinates of the mesh vertices.

minval, maxval - The minimum and maximum values on the color scale. A value of minval will appear completely blue and a value of maxval will appear completely red. If no values are specified, minval will default to the lowest value of the variable over the area to be plotted and maxval will default to the highest value of the variable over the area to be plotted.

resolution - The distance between the centers of adjacent plot cells. Lower values will result in higher quality images. If no value is specified, the resolution will be set such that there will be approximately 10,000 plot cells.

filename - The name of the file to store the resulting image in. An extension is not needed here - the .eps extension will automatically be added. If no filename is specified, the image will be stored in the file `temp.eps` so it can be viewed.

viewcmd - The OS command used to access your graphics viewer. If specified, the graphics viewer will be automatically opened
to view your plot. Default is 'gv'.

xlabel - The label to use for the X axis. Default is 'X values'.

ylabel - The label to use for the Y axis. Default is 'Y values'.

valuelabel - The label to use for the value of the variable on the scale. Default is the variable's name.

gridcorrect - The amount, in multiples of the resolution, that adjacent plot cells overlap. This overlap eliminates the 'grid effect'. If you still
see a grid of white space, increase this value. Default value is 0.3.

graphheight - The physical height of the graph. Default is 12.

graphwidth - The physical width of the graph. Default is 12.

mask - Used to 'blank out' specified areas of the graph, useful if you have a mesh whose boundaries do not form a rectangle. The format for mask is an N-by-4 array, where N is the number of different rectangles you want to blank out. Each row represents a separate blanked-out rectangle: the first two numbers being the X and Y coordinates of the lower left hand corner of said rectangle, and the last two being the X and Y coordinates of the upper right hand corner. For example, mask = [[0, 0, 3, 2], [7, 8, 10, 10]] indicates that the area inside the rectangle with lower left corner (0, 0) and upper left corner (3, 2) should be blanked out, as should the area inside the rectangle with lower left coordinate (7, 8) and upper left coordinate (10, 10).

The following are keyword arguments in the __init__() method: (Note that these work for PyxViewer only, not Grid2DPyxViewer)

showpercent - If set to 1, this displays the percentage complete at 5% intervals as the plot progresses. Default value is 1.

showtime - If set to 1, this displays the approximate time it will take to plot it. Default value is 1.

Test cases:

   >>> from fipy.meshes.grid2D import Grid2D
   >>> import fipy.variables.cellVariable
   >>> mesh = Grid2D(1.0, 1.0, 2, 2)
   >>> myvar = fipy.variables.cellVariable.CellVariable(mesh)
   >>> myvar[0] = 1.0
   >>> myvar[1] = 2.0
   >>> myvar[2] = 3.0
   >>> myvar[3] = 4.0
   >>> myviewer = Grid2DPyxViewer(myvar, showpercent = 0, showtime = 0)
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
   >>> myviewer = Grid2DPyxViewer(myvar, showpercent = 0, showtime = 0)
   >>> testlist = myviewer.plot(minx = 0.0, maxx = 2.0, miny = 0.0, maxy = 2.0, resolution = 1.0, minval = 1.0, maxval = 5.0, gridcorrect = 0.0, returnlist = 1, viewcmd = None) 
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
   >>> testlist = myviewer.plot(minx = 0.0, maxx = 2.0, miny = 0.0, maxy = 2.0, resolution = 1.0, gridcorrect = 0.0, returnlist = 1, viewcmd = None) 
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
   >>> myviewer = PyxViewer(myvar, showpercent = 0, showtime = 0)
   >>> array = myviewer.getValueMatrix(0.0, 2.0, 0.0, 2.0, 0.5)
   >>> testlist = array.tolist()
   >>> print testlist
   [[1.0, 1.0, 3.0, 3.0], [1.0, 1.0, 3.0, 3.0], [2.0, 2.0, 4.0, 4.0], [2.0, 2.0, 4.0, 4.0]]

   >>> from fipy.meshes.numMesh.grid3D import Grid3D
   >>> import fipy.variables.cellVariable
   >>> mesh = Grid3D(1.0, 1.0, 1.0, 2, 2, 2)
   >>> myvar = fipy.variables.cellVariable.CellVariable(mesh)
   >>> myvar[0] = 1.0
   >>> myvar[1] = 2.0
   >>> myvar[2] = 3.0
   >>> myvar[3] = 4.0
   >>> myvar[4] = 5.0
   >>> myvar[5] = 6.0
   >>> myvar[6] = 7.0
   >>> myvar[7] = 8.0
   >>> myviewer = Grid3DPyxViewer(myvar, zvalue = 1.0)
   >>> array = myviewer.getValueMatrix(0.0, 2.0, 0.0, 2.0, 0.5)
   >>> testlist = array.tolist()
   >>> print testlist
   [[2.25, 3.25, 4.25, 5.25], [2.75, 3.75, 4.75, 5.75], [3.25, 4.25, 5.25, 6.25], [3.75, 4.75, 5.75, 6.75]]
   
   >>> from fipy.meshes.numMesh.grid3D import Grid3D
   >>> import fipy.variables.cellVariable
   >>> mesh = Grid3D(1.0, 1.0, 1.0, 2, 2, 1)
   >>> myvar = fipy.variables.cellVariable.CellVariable(mesh)
   >>> myvar[0] = 1.0
   >>> myvar[1] = 2.0
   >>> myvar[2] = 3.0
   >>> myvar[3] = 4.0
   >>> myviewer = Grid3DPyxViewer(myvar, zvalue = 0.75)
   >>> array = myviewer.getValueMatrix(0.0, 2.0, 0.0, 2.0, 0.5)
   >>> testlist = array.tolist()
   >>> print testlist
   [[0.25, 1.25, 2.25, 3.25], [0.75, 1.75, 2.75, 3.75], [1.25, 2.25, 3.25, 4.25], [1.75, 2.75, 3.75, 4.75]]
   
"""

__docformat__ = 'restructuredtext'

import Numeric
import os
import time

if not os.environ.has_key('FIPY_NOPYX'):
    import pyx
    import pyx.color

    pyx.text.set(fontmaps="psfonts.cmz")

class PyxViewer:
    
    def __init__(self, variable, showpercent = 1, showtime = 1):
        self.var = variable
        self.showpercent = showpercent
        self.showtime = showtime

## ------------------------------------------------------------------------------------

    def plot(self, returnlist = 0, debug=0, minx=None, maxx=None, miny=None, maxy=None, minval=None, maxval=None, resolution = None, filename = None, viewcmd = "gv", xlabel = "X values", ylabel = "Y values", valuelabel = None, mask = None, gridcorrect = 0.3, graphheight = 12, graphwidth = 12):
        
        ## initialize variables
        
        starttime = time.time()
        thepalette = pyx.color.palette.ReverseRainbow
        theMesh = self.var.getMesh()
        vertexCoords = theMesh.getVertexCoords()
        
        if(minx == None):
            minx = vertexCoords[Numeric.argmin(vertexCoords, axis = 0)[0], 0]
        if(miny == None):
            miny = vertexCoords[Numeric.argmin(vertexCoords, axis = 0)[1], 1]
        if(maxx == None):
            maxx = vertexCoords[Numeric.argmax(vertexCoords, axis = 0)[0], 0]
        if(maxy == None):
            maxy = vertexCoords[Numeric.argmax(vertexCoords, axis = 0)[1], 1]
            
        ## calculate the resolution. If there is no resolution set, calculate a "default" resolution that will lead to 1,000 points being plotted.
        ## If a resolution is given, use that.
        
        if(resolution == None):
            resolution = ((((maxx - minx)*(maxy - miny))/10000.0) ** 0.5)

        ## calculate the number of points to plot in the X direction (xsize) and the Y direction (ysize). The total number of points plotted
        ## will be equal to (xsize * ysize).
        xsize = ((maxx - minx) / resolution) - 0.5
        ysize = ((maxy - miny) / resolution) - 0.5
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
        if(debug == 1):
            print xarr
            print yarr
            print valarr
        ## This section takes the values above and generates the list to put into the pyx.graph.plot method.
        ## If maxval and minval are not given, generate them.
        xlist = Numeric.reshape(xarr, (Numeric.size(xarr),))
        ylist = Numeric.reshape(yarr, (Numeric.size(yarr),))
        vallist = Numeric.reshape(valarr, (Numeric.size(valarr),))
        if(minval == None):
            minval = min(vallist)
        if(maxval == None):
            maxval = max(vallist)
        if(debug == 1):
            print vallist
            print minval
            print maxval
        minval = float(minval)
        maxval = float(maxval)
        vallist = vallist.astype(Numeric.Float)
	if maxval == minval:
	    vallist = vallist - minval
	else:
	    vallist = (vallist - minval) / (maxval - minval)
        minxlist = xlist - (0.5 * resolution)
        maxxlist = xlist + ((0.5 + gridcorrect) * resolution)
        minylist = ylist - (0.5 * resolution)
        maxylist = ylist + ((0.5 + gridcorrect) * resolution)
        baselist = Numeric.array([minxlist, maxxlist, minylist, maxylist, vallist])
        displaylist = Numeric.transpose(baselist)
        ## If the user did not specify a filename but does want to view the picture, create a file anyway so the picture can be viewed.
        if(filename == None):
            if(viewcmd != None):
                filename = "temp"
        displayinput = displaylist.tolist()
        ## Create the graph and plot it.
        if(filename != None):
            mycanvas = pyx.canvas.canvas()
            mygraph = pyx.graph.graphxy(height = graphheight, width = graphwidth, xpos = 0, ypos = 0, x = pyx.graph.axis.linear(min = minx, max = maxx, title = xlabel), y = pyx.graph.axis.linear(min = miny, max = maxy, title = ylabel))
            mygraph.plot(pyx.graph.data.list(displayinput, xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5), pyx.graph.style.rect(thepalette))
            mycanvas.stroke(mygraph)
            if(valuelabel == None):
                valuelabel = self.var.name
            scalegraph = pyx.graph.graphxy(height = 12, width = 1, xpos = graphwidth + 2, ypos = 0, y = pyx.graph.axis.linear(min = minval, max = maxval, title = valuelabel))
            scaleminx = Numeric.zeros((100,))
            scalemaxx = Numeric.ones((100,))
            scaleminy = Numeric.fromfunction(lambda num: minval + (maxval - minval)*(num / 100.0), (100,))
            scalemaxy = Numeric.fromfunction(lambda num: minval + (maxval - minval)*((num + (1 + gridcorrect)) / 100.0), (100,))
            scalecolors = Numeric.fromfunction(lambda num: (num + 0.5) / 100.0, (100,))
            scalearray = Numeric.array([scaleminx, scalemaxx, scaleminy, scalemaxy, scalecolors])
            scalearray = Numeric.transpose(scalearray)
            scalelist = scalearray.tolist()
            scalegraph.plot(pyx.graph.data.list(scalelist, xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5), pyx.graph.style.rect(thepalette))
            mycanvas.stroke(scalegraph)
            ## do the mask
            if(mask != None):
                ## the following two functions convert "grid" coordinates (the coordinate system used by the variables and mesh) to "display" coordinates (the coordinates used by the canvas)
                xgridtodisplay = lambda x: min(((x - minx)/(maxx - minx)) * graphwidth, graphwidth) 
                ygridtodisplay = lambda y: min(((y - miny)/(maxy - miny)) * graphheight, graphheight)
                for element in mask:
                    lowx = xgridtodisplay(element[0])
                    lowy = ygridtodisplay(element[1])
                    highx = xgridtodisplay(element[2])
                    highy = ygridtodisplay(element[3])
                    blankpart = pyx.path.rect(lowx, lowy, (highx - lowx), (highy - lowy))
                    mycanvas.stroke(blankpart, [pyx.deco.filled([pyx.color.grey(0.95)])])
            mycanvas.writeEPSfile(filename)
            if(viewcmd != None):
                os.system(viewcmd + " " + filename + ".eps &")
        if(returnlist == 1):
            return displayinput
        else:
            return None
## ---------------------------------------------------------------------------------------------------
    def generateArray(self):
        """ generates the array of variable values
        """
        self.theArray = Numeric.array(self.var)
## ---------------------------------------------------------------------------------------------------
            
    def getValueMatrix(self, minx, maxx, miny, maxy, resolution):
        self.generateArray() 
        xsize = ((maxx - minx) / resolution) - 0.5
        ysize = ((maxy - miny) / resolution) - 0.5
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
        xsize = x.shape[0]
        ysize = x.shape[1]
        totsize = xsize * ysize
        if(totsize > 10000):
            estamt = 100
        else:
            estamt = (totsize / 100)
        if estamt == 0:
            estamt = 1
        resArray = Numeric.zeros((totsize,))
        resArray = resArray.astype(Numeric.Float)
        xlist = Numeric.reshape(x, (totsize,))
        ylist = Numeric.reshape(y, (totsize,))
        cellCenters = self.var.getMesh().getCellCenters()
        cellXvalues = Numeric.transpose(cellCenters)[0]
        cellYvalues = Numeric.transpose(cellCenters)[1]
        meshsize = Numeric.size(cellXvalues)
        plotCellsPerLoop = 1
        cellXvalues = Numeric.reshape(cellXvalues, (meshsize, 1))
        cellYvalues = Numeric.reshape(cellYvalues, (meshsize, 1))
        currPlot = 0
        lastDisplayedPercent = 0
        displayInterval = 5
        if(self.showpercent == 1):
            print "Entering process now."
        starttime = time.time()
        while (currPlot + plotCellsPerLoop < totsize):
            Xoffsets = cellXvalues - xlist[currPlot : currPlot + plotCellsPerLoop]
            Yoffsets = cellYvalues - ylist[currPlot : currPlot + plotCellsPerLoop]
            a = Xoffsets ** 2
            b = Yoffsets ** 2
            squaredDists = a + b
            cellIDarray = Numeric.argmin(squaredDists, axis = 0)
            resArray[currPlot: currPlot + plotCellsPerLoop] = Numeric.take(self.theArray, cellIDarray)
            currPlot = currPlot + plotCellsPerLoop
            if(self.showtime == 1):
                if(currPlot == estamt):
                    estTime = ((time.time() - starttime) * totsize) / estamt
                    estTime = int(estTime)
                    print "Estimated time to completion:", estTime, "seconds." 
            if(self.showpercent == 1):            
                percentComplete = (currPlot * 100) / totsize
                if(percentComplete - lastDisplayedPercent > displayInterval):
                    lastDisplayedPercent = lastDisplayedPercent + displayInterval
                    print lastDisplayedPercent, "percent complete."
        Xoffsets = cellXvalues - xlist[currPlot:]
        Yoffsets = cellYvalues - ylist[currPlot:]
        a = Xoffsets ** 2
        b = Yoffsets ** 2
        squaredDists = a + b
        cellIDarray = Numeric.argmin(squaredDists, axis = 0)
        resArray[currPlot:] = Numeric.take(self.theArray, cellIDarray)
        if(self.showpercent == 1):
            print "100 percent complete."
        resArray = Numeric.reshape(resArray, (xsize, ysize))
        return resArray

##------------------##------------------##------------------------------## ----------------##----------------##----------

InvalidClassException = "Wrong class"

class Grid2DPyxViewer(PyxViewer):
    def __init__(self, variable, showpercent = 0, showtime = 0):
        if(variable.getMesh().__class__.__name__ != "Grid2D"):
            raise InvalidClassException, "Grid2DPyxViewer can be used only with Grid2D meshes"
        else:
            self.var = variable
    
    def getValue(self, x, y):
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
        if(width == 1):
            xnumcells = Numeric.zeros(xnumcells.shape)
            Tvalues = Numeric.zeros(xnumcells.shape)
        if(height == 1):
            ynumcells = Numeric.zeros(ynumcells.shape)
            Uvalues = Numeric.zeros(ynumcells.shape)
        Y1coeffs = (1 - Tvalues)*(1 - Uvalues)
        Y2coeffs = Tvalues * (1 - Uvalues)
        Y3coeffs = Tvalues * Uvalues
        Y4coeffs = (1 - Tvalues) * Uvalues
        cellIDarray = (ynumcells * width) + xnumcells
        resArray = Numeric.take(self.theArray, cellIDarray) * Y1coeffs   ## contribution from lower left point
        if(width != 1):
            resArray = resArray + Numeric.take(self.theArray, cellIDarray + 1) * Y2coeffs   ## contribution from lower right point
        if(height != 1):
            resArray = resArray + Numeric.take(self.theArray, cellIDarray + (width + 1)) * Y3coeffs   ## contribution from upper right point
        if(width != 1 and height != 1):
            resArray = resArray + Numeric.take(self.theArray, cellIDarray + width) * Y4coeffs   ## contribution from upper left point
        return resArray


##----------------------##---------------------------##------------------------##-----------------------##---------------------##

ValueOutOfRangeException = "Value out of range"

class Grid3DPyxViewer(Grid2DPyxViewer):
    def __init__(self, variable, showpercent = 0, showtime = 0, zvalue = 0):
        if(variable.getMesh().__class__.__name__ != "Grid3D"):
            raise InvalidClassException, "Grid3DPyxViewer can be used only with Grid3D meshes"
        else:
            self.var = variable
            if(zvalue < 0) or(zvalue > (self.var.getMesh().__getstate__()["nz"] * self.var.getMesh().__getstate__()["dz"])):
                raise ValueOutOfRangeException, "Z value is out of the range of the mesh"
            self.zvalue = zvalue

    def generateArray(self):
        totarray = Numeric.array(self.var)
        if(self.var.getMesh().__getstate__()["nz"] == 1):
            self.theArray = totarray
            return None
        zcells = (self.zvalue / self.var.getMesh().__getstate__()["dz"]) - 0.5
        znumcells = int(zcells)
        if(znumcells == (self.var.getMesh().__getstate__()["nz"] - 1)):
            znumcells = znumcells - 1
        resid = zcells - znumcells
        totarray = Numeric.array(self.var)
        persheet = (self.var.getMesh().__getstate__()["nx"]) * (self.var.getMesh().__getstate__()["ny"])
        bottom = totarray[persheet * znumcells : persheet * (znumcells + 1)]
        top = totarray[persheet * (znumcells + 1) : persheet * (znumcells + 2)]
        result = (top * resid) + (bottom * (1.0 - resid))
        self.theArray = result

##--------------------------##---------------------##----------------------------##----------------------##-------------------------##
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
