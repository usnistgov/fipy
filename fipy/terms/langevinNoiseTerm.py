#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "langevinNoiseTerm.py"
 #                                    created: 7/26/05 {8:35:17 AM} 
 #                                last update: 7/26/05 {12:24:27 PM} 
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
 # protection and is in the public domain.  langevinNoiseTerm.py
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
 #  2005-07-26 JEG 1.0 original
 # ###################################################################
 ##

"""
	>>> from fipy.meshes.grid2D import Grid2D
        >>> v = LangevinNoiseTerm(mesh = Grid2D(nx = 100, ny = 100), stddev = 2)
        >>> from fipy import viewers
        >>> viewer = viewers.make(vars = v, limits = {'datamin': -3, 'datamax':3})
        
        >>> from fipy.tools.numerix import *
        
        >>> def histogram(a, bins):
        ...     n = searchsorted(sort(a), bins)
        ...     n = concatenate([n, [len(a)]])
        ...     return n[1:] - n[:-1]
        
        >>> from fipy.meshes.grid1D import Grid1D
        >>> from fipy.variables.cellVariable import CellVariable
        >>> hist = CellVariable(mesh = Grid1D(nx = 200, dx = .1) + (-10.,))
        >>> histoplot = viewers.make(vars = hist)
        
        >>> for i in range(1000):
        ...     viewer.plot()
        ...     hist[:] = histogram(v(), bins = hist.getMesh().getCellCenters()[:,0])
        ...     histoplot.plot()

"""

__docformat__ = 'restructuredtext'

from RandomArray import normal

from fipy.variables.cellVariable import CellVariable

class LangevinNoiseTerm(CellVariable):
    def __init__(self, mesh, name = '', mean = 0., stddev = 1., unit = None):
        self.mean = mean
        self.stddev = stddev
        CellVariable.__init__(self, mesh = mesh, name = name, 
                              value = 0, unit = unit, hasOld = 0)
    
    def _refresh(self):
        """
        We *always* want the random variable to be stale, so that we
        get new random values whenever we look at it.
        """
        CellVariable._refresh(self)
        self._markStale()
        
    def _calcValue(self):
        self.value = normal(self.mean, self.stddev, 
                            shape = [self.getMesh().getNumberOfCells()])

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 

