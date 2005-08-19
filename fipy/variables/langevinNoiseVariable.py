#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "langevinNoiseVariable.py"
 #                                    created: 7/26/05 {8:35:17 AM} 
 #                                last update: 8/19/05 {12:46:36 PM} 
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
        >>> mean = 0.
        >>> stddev = 2.
        
	>>> from fipy.meshes.grid2D import Grid2D
        >>> v = LangevinNoiseTerm(mesh = Grid2D(nx = 100, ny = 100), mean = mean, stddev = stddev)
        >>> vdg = v.getFaceGrad().getDivergence()
        >>> from fipy import viewers
        >>> viewer = viewers.make(vars = v)
        
        >>> from fipy.tools.numerix import *
        
        >>> def histogram(a, bins):
        ...     n = searchsorted(sort(a), bins)
        ...     n = concatenate([n, [len(a)]])
        ...     dx = bins[1:] - bins[:-1]
        ...     return (n[1:] - n[:-1]) / concatenate([dx, [dx[-1]]]) / float(len(a))
        
        >>> from fipy.meshes.grid1D import Grid1D
        >>> from fipy.variables.cellVariable import CellVariable
        >>> histomesh = Grid1D(nx = 600, dx = .1) + (-30.,)
        >>> hist = CellVariable(mesh = histomesh)
        >>> gauss = CellVariable(mesh = histomesh)
        >>> x = histomesh.getCellCenters()[...,0]
        >>> gauss.setValue((1/(stddev * sqrt(2*pi))) * exp(-(x - mean)**2 / (2 * stddev**2)))
        >>> histoplot = viewers.make(vars = (hist, gauss))
        
        >>> for i in range(10):
        ...     v.updateOld()
        ...     print vdg.getCellVolumeAverage()
        ...     0.
        ...     hist[:] = histogram(v(), bins = x)
        ...     if __name__ == '__main__':
        ...         viewer.plot()
        ...         histoplot.plot()

"""

__docformat__ = 'restructuredtext'

from RandomArray import normal

from fipy.variables.cellVariable import CellVariable

class LangevinNoiseVariable(CellVariable):
    def __init__(self, mesh, name = '', mean = 0., stddev = 1., unit = None, hasOld = 1):
        self.mean = mean
        self.stddev = stddev
        CellVariable.__init__(self, mesh = mesh, name = name, 
                              value = 0, unit = unit, hasOld = hasOld)
        if hasOld:
            self._markStale()
    
    def copy(self):
        cpy = self.__class__(
            mesh = self.mesh, 
            name = self.name + "_old", 
            mean = self.mean,
            stddev = self.stddev,
            unit = self.getUnit(),
            hasOld = 0)
        cpy.setValue(self.getValue())
        return cpy

    def updateOld(self):
        CellVariable.updateOld(self)
        self._markStale()
        
    def _calcValue(self):
        self.value = normal(self.mean, self.stddev, 
                            shape = [self.getMesh().getNumberOfCells()])

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 

