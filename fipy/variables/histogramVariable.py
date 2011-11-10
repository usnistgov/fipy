## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "histogramVariable.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  histogramVariable.py
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
 #  See the file "license.terms" for information on usage and  redistribution of
 #  this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable
from fipy.meshes import Grid1D
from fipy.tools import numerix

__all__ = ["HistogramVariable"]

class HistogramVariable(CellVariable):
    def __init__(self, distribution, dx = 1., nx = None, offset = 0.):
        r"""
        Produces a histogram of the values of the supplied distribution.
        
        :Parameters:
            
            - `distribution`: The collection of values to sample.
            - `dx`: the bin size
            - `nx`: the number of bins
            - `offset`: the position of the first bin
        """
        CellVariable.__init__(self, mesh = Grid1D(dx = dx, nx = nx) + (offset,))
        self.distribution = self._requires(distribution)
        
    def _calcValue(self):
        l = len(self.distribution)
        bins = self.mesh.cellCenters[0]
        n = numerix.searchsorted(numerix.sort(self.distribution), bins)
        n = numerix.concatenate([n, [l]])
        dx = bins[1:] - bins[:-1]
        return (n[1:] - n[:-1]) / numerix.concatenate([dx, [dx[-1]]]) / float(l)
