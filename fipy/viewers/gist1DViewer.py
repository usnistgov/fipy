#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "gist1DViewer.py"
 #                                    created: 11/10/03 {2:48:25 PM} 
 #                                last update: 1/20/04 {4:16:40 PM} 
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
 
import gist

from fivol.viewers.gistViewer import GistViewer

class Gist1DViewer(GistViewer):
    
    def __init__(self, vars = None, title = None, minVal=None, maxVal=None, xlog = 0, ylog = 0):
        self.vars = list(vars)
	self.xlog = xlog
	self.ylog = ylog
	if title is None and len(self.vars) == 1:
	    title = self.vars[0].name
	else:
	    title = ''
        GistViewer.__init__(self, minVal, maxVal, title = title)

    def getArrays(self):
	arrays = ()
	for var in self.vars:
	    nx,ny = var.getMesh().getShape()
## 	    arrays += Numeric.reshape(var.getNumericValue(),(ny,nx)),
	    arrays += var.getNumericValue(),
	return arrays

    def plot(self, minVal=None, maxVal=None):
	gist.window(self.id, wait= 1)
	gist.pltitle(self.title)
	gist.animate(1)
	
## 	if minVal is None:
## 	    if self.minVal is None:
## 		minVal = Numeric.minimum.reduce(array.flat)
## 	    else:
## 		minVal = self.minVal
## 	if maxVal is None:
## 	    if self.maxVal is None:
## 		maxVal = Numeric.minimum.reduce(array.flat)
## 	    else:
## 		maxVal = self.maxVal
## 
## ## 	min = Numeric.minimum.reduce(array.flat)
## ## 	max = Numeric.minimum.reduce(array.flat)
## 	if maxVal == minVal:
## 	    maxVal = minVal + 0.01

	for array in self.getArrays():
	    gist.plg(array)
	    gist.logxy(self.xlog, self.ylog)
	gist.fma()
