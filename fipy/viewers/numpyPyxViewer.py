
#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "pyxviewer.py"
 #                                    created: 6/25/04 {3:17:21 PM} 
 #                                last update: 11/3/04 {11:12:18 AM} 
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

`NumpyPyxViewer` - A simple PyxViewer for 2D arrays

"""

import Numeric

class NumpyPyxViewer:
    def __init__(self, array, dx = 1., dy = 1., plotWidth = 10, viewCommand = 'gv', maxVal = 1., minVal = 0.):
        self.array = Numeric.array(array)
        self.dx = dx
        self.dy = dy
        (self.nx, self.ny) = self.array.shape
        self.plotWidth = plotWidth
        self.viewCommand = viewCommand
        self.maxVal = maxVal
        self.minVal = minVal

        self.data = Numeric.zeros((self.nx * self.ny, 5), 'd')
        indices = Numeric.indices(self.array.shape)
        self.data[:,0] = indices[1].flat * self.dx
        self.data[:,1] = indices[1].flat * self.dx + self.dx
        self.data[:,2] = indices[0].flat * self.dy
        self.data[:,3] = indices[0].flat * self.dy + self.dy
        
    def plot(self, fileName = 'tmp'):
        Lx = self.dx * self.nx
        Ly = self.dy * self.ny

	import pyx
	
        g = pyx.graph.graphxy(height = self.plotWidth * Ly / Lx, width = self.plotWidth,
                              x = pyx.graph.axis.linear(min = 0, max = Lx),
                              y = pyx.graph.axis.linear(min = 0, max = Ly))

        clippedarray = Numeric.where(self.array > self.maxVal, self.maxVal, self.array)
        clippedarray = Numeric.where(clippedarray < self.minVal, self.minVal, clippedarray) 

        clippedarray = (clippedarray - self.minVal) / (self.maxVal - self.minVal)
    
        self.data[:,4] = clippedarray.flat

        g.plot(pyx.graph.data.list(self.data.tolist(), xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5),
               pyx.graph.style.rect(pyx.color.palette.Rainbow))

        g.dodata()

        g.writeEPSfile(fileName)

	import os
        os.system(self.viewCommand + " " + fileName + ".eps &")


