## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "exponentialNoiseVariable.py"
 #                                    created: 8/27/05 {9:26:58 AM} 
 #                                last update: 8/27/05 {9:53:44 AM} 
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
 # protection and is in the public domain.  exponentialNoiseVariable.py
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

from RandomArray import exponential

from fipy.variables.noiseVariable import NoiseVariable

class ExponentialNoiseVariable(NoiseVariable):
    r"""
    Represents an exponential distribution of random numbers with the probability
    distribution
    
    .. raw:: latex
    
       \[ \mu^{-1} e^{-\frac{x}{\mu}} \]
       
       with a mean parameter $\mu$.

    We generate noise on a uniform cartesian mesh
           
           >>> from fipy.variables.variable import Variable
           >>> mean = Variable()
           >>> from fipy.meshes.grid2D import Grid2D
           >>> noise = ExponentialNoiseVariable(mesh = Grid2D(nx = 100, ny = 100), mean = mean)
           
    We histogram the root-volume-weighted noise distribution
    
           >>> from fipy.variables.histogramVariable import HistogramVariable
           >>> histogram = HistogramVariable(distribution = noise, dx = 0.1, nx = 100)
           
    and compare to a Gaussian distribution
    
           >>> from fipy.variables.cellVariable import CellVariable
           >>> expdist = CellVariable(mesh = histogram.getMesh())
           >>> x = histogram.getMesh().getCellCenters()[...,0]
           
           >>> if __name__ == '__main__':
           ...     from fipy import viewers
           ...     viewer = viewers.make(vars = noise, limits = {'datamin': 0, 'datamax': 5})
           ...     histoplot = viewers.make(vars = (histogram, expdist), limits = {'datamin': 0, 'datamax': 1.5})
           
           >>> from fipy.tools.numerix import arange, exp
           
           >>> mean.setValue(1.5)
           >>> viewer.plot(filename = "exp.cgm")
           >>> expdist.setValue((1/mean)*exp(-x/mean))
           >>> from fipy.viewers.tsvViewer import TSVViewer
           >>> TSVViewer(vars = (histogram, expdist)).plot(filename = "exp.tsv")
           
           >>> for mu in arange(0.5,3,0.5):
           ...     mean.setValue(mu)
           ...     import sys
           ...     print >>sys.stderr, "mean: %g" % mean
           ...     expdist.setValue((1/mean)*exp(-x/mean))
           ...     if __name__ == '__main__':
           ...         viewer.plot()
           ...         histoplot.plot()

           >>> print abs(noise.getFaceGrad().getDivergence().getCellVolumeAverage()) < 5e-15
           1

    Note that the noise exhibits larger amplitude in the small cells than in the large ones

    .. image:: fipy/variables/exp.jpg
      :scale: 25
      :align: center

    but that the root-volume-weighted histogram is Gaussian.

    .. image:: fipy/variables/exp-histogram.pdf
      :scale: 25
      :align: center

    """
    def __init__(self, mesh, mean, name = '', hasOld = 0):
        r"""
        :Parameters:
            - `mesh`: The mesh on which to define the noise.
            - `mean`: The mean of the distribution
            
              .. raw:: latex
              
                 $\mu$.
        """
        NoiseVariable.__init__(self, mesh = mesh, name = name, hasOld = hasOld)
        self.mean = self._requires(mean)
    
    def _calcValue(self):
        self.value = exponential(mean = self.mean, 
                                 shape = [self.getMesh().getNumberOfCells()])

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 



