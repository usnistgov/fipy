#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gistViewer.py"
 #                                    created: 11/10/03 {2:48:25 PM} 
 #                                last update: 3/4/05 {4:28:18 PM} { 2:45:36 PM}
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

import os

import Numeric
import gist

from fipy.viewers.viewer import Viewer

class GistViewer(Viewer):
    
    id=0
    
    def __init__(self, vars, limits = None, title = None, dpi = 75):
        """
        :Parameters:
          - `vars`: a `Variable` or tuple of `Variable` objects to plot
          - `limits`: a dictionary with possible keys `xmin`, `xmax`, 
                      `ymin`, `ymax`, `zmin`, `zmax`, `datamin`, `datamax`.
                      A 1D Viewer will only use `xmin` and `xmax`, a 2D viewer 
                      will also use `ymin` and `ymax`, and so on. 
                      All viewers will use `datamin` and `datamax`. 
                      Any limit set to a (default) value of `None` will autoscale.
          - `title`: displayed at the top of the Viewer window
          - `dpi`: the dot-per-inch resolution of the display
        """
        Viewer.__init__(self, vars = vars, limits = limits, title = title)
        
        self.mesh = self.vars[0].getMesh()

        self.id = GistViewer.id 
	GistViewer.id += 1
        
        gist.window(self.id, wait = 1, dpi = dpi, display = '')

    def getLimit(self, key):
        limit = Viewer.getLimit(self, key = key)
        if limit is None:
            limit = 'e'
            
        return limit
        
    def plot(self):
        pass


        
