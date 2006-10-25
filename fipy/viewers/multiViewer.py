## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "multiViewer.py"
 #                                    created: 12/17/05 {6:28:03 PM} 
 #                                last update: 10/25/06 {4:10:46 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  multiViewer.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
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
 #  2005-12-17 JEG 1.0 original
 # ###################################################################
 ##

from fipy.viewers.viewer import Viewer

class MultiViewer(Viewer):
    """
    Treat a collection of different viewers (such for different 2D plots 
    or 1D plots with different axes) as a single viewer that will `plot()` 
    all subviewers simultaneously.
    """
    def __init__(self, viewers):
        if type(viewers) not in [type([]), type(())]:
            viewers = [viewers]
        self.viewers = viewers

    def setLimits(self, limits):
        for viewer in self.getViewers():
            viewer.setLimits(limits)
            
    def plot(self):
        for viewer in self.getViewers():
            viewer.plot()
            
    def getViewers(self):
        return self.viewers
