#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gist1DViewer.py"
 #                                    created: 11/10/03 {2:48:25 PM} 
 #                                last update: 11/2/04 {7:38:00 PM} 
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

import Numeric
 
import os
try:
    import gist
except:
    print "Unable to load gist"

from fipy.viewers.gistViewer import GistViewer

class Gist1DViewer(GistViewer):
    
    def __init__(self, vars = None, title = None, limits = None, xlog = 0, ylog = 0, style = "work.gs"):
        self.vars = list(vars)
	self.xlog = xlog
	self.ylog = ylog
	self.style = style
	if title is None and len(self.vars) == 1:
	    title = self.vars[0].name
	else:
	    title = ''
        GistViewer.__init__(self, limits = limits, title = title)

    def getArrays(self):
	arrays = ()
	for var in self.vars:
##	    arrays += var.getNumericValue(),
            arrays += (Numeric.array(var),)
	return arrays
	
    def plotArrays(self):
## 	gist.plsys(0)
	for array in self.getArrays():
            if len(array.shape) > 1:
                gist.plg(array[:,1], array[:,0])
            else:
                gist.plg(array)
	gist.logxy(self.xlog, self.ylog)

    def plot(self, minVal=None, maxVal=None):
	gist.window(self.id, wait= 1, style = self.style)
	gist.pltitle(self.title)
	gist.animate(1)

	if self.limits != None:
	    gist.limits(self.limits[0], self.limits[1], self.limits[2], self.limits[3])
	    
	self.plotArrays()
	    
	gist.fma()
