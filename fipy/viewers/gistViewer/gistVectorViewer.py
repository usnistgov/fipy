#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gistViewer.py"
 #                                    created: 11/10/03 {2:48:25 PM} 
 #                                last update: 7/6/05 {4:39:36 PM} { 2:45:36 PM}
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

from fipy.viewers.gistViewer import GistViewer

class _GistVectorViewer(GistViewer):
    
    def __init__(self, var = None, title = '', grid = 1):
	self.var = var
	GistViewer.__init__(self, title, grid)
	
    def plot(self, filename = None):
	import gist

        gist.window(self.id, wait = 1)
	gist.pltitle(self.title)
        gist.animate(1)
        gist.palette(self.palette)
## 	gist.gridxy(self.grid)
	
	centers = self.var.getMesh().getFaceCenters()
	
	gist.plmesh(Numeric.array([centers[:,1],centers[:,1]]), Numeric.array([centers[:,0],centers[:,0]]))

	print self.var
	
	vx = Numeric.array(self.var[:,0])
	vy = Numeric.array(self.var[:,1])
	
        gist.plv(Numeric.array([vy,vy]), Numeric.array([vx,vx]),scale=0.01)

        if filename is not None:
            
            gist.hcp_file(filename)
            gist.hcp()

        gist.fma()
        
    def getArray(self):
        pass
        
