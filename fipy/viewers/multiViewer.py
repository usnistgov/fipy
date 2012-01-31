## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "multiViewer.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
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
 # ###################################################################
 ##

from fipy.viewers.viewer import AbstractViewer
from fipy.tools.decorators import getsetDeprecated

__all__ = ["MultiViewer"]

class MultiViewer(AbstractViewer):
    """
    Treat a collection of different viewers (such for different 2D plots 
    or 1D plots with different axes) as a single viewer that will `plot()` 
    all subviewers simultaneously.
    """
    def __init__(self, viewers):
        """
        :Parameters:
          viewers : list
            the viewers to bind together
        """
        if type(viewers) not in [type([]), type(())]:
            viewers = [viewers]
        self.viewers = viewers

    def setLimits(self, limits={}, **kwlimits):
        kwlimits.update(limits)
        for viewer in self.viewers:
            viewer.setLimits(**kwlimits)
            
    def plot(self):
        for viewer in self.viewers:
            viewer.plot()
            
    @getsetDeprecated
    def getViewers(self):
        return self.viewers
