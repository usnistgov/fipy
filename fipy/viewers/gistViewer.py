#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gistViewer.py"
 #                                    created: 11/10/03 {2:48:25 PM} 
 #                                last update: 9/1/04 {6:38:32 PM} { 2:45:36 PM}
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
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

import Numeric

import os
if not os.environ.has_key('FIPY_NOGIST'):
    import gist
    import colorbar

class GistViewer:
    
    id=0
    
    def __init__(self, minVal = None, maxVal = None, title = '', palette = 'heat.gp', grid = 1):
	self.minVal = minVal
        self.maxVal = maxVal
	self.title = title
        self.id = GistViewer.id 
	GistViewer.id += 1
        self.palette = palette
        self.grid = grid
        
    def plot(self, minVal = None, maxVal = None, fileName = None):
	array = self.getArray()

        gist.window(self.id, wait = 1)
	gist.pltitle(self.title)
        gist.animate(1)
        gist.palette(self.palette)
	gist.gridxy(self.grid)
	
	if minVal is None:
	    if self.minVal is None:
		minVal = Numeric.minimum.reduce(array.flat)
	    else:
		minVal = self.minVal
                
	if maxVal is None:
	    if self.maxVal is None:
		maxVal = Numeric.minimum.reduce(array.flat)
	    else:
		maxVal = self.maxVal
        
	if maxVal == minVal:
	    maxVal = minVal + 0.01

        gist.pli(array, cmin = minVal, cmax = maxVal)
	colorbar.color_bar(minz = minVal, maxz = maxVal, ncol=240, zlabel = 'fred')

        if fileName is not None:
            
            gist.hcp_file(fileName)
            gist.hcp()

        gist.fma()
        
    def getArray(self):
        pass
        
