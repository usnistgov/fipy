## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "gammaNoiseVariable.py"
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
 # protection and is in the public domain.  gammaNoiseVariable.py
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

from fipy.tools.numerix import random

from fipy.variables.noiseVariable import NoiseVariable

class GammaNoiseVariable(NoiseVariable):
    r"""
    Represents a gamma distribution of random numbers with the probability
    distribution
    
    .. math::
    
       x^{\alpha - 1}\frac{\beta^\alpha e^{-\beta x}}{\Gamma(\alpha)}
       
    with a shape parameter :math:`\alpha`, a rate parameter :math:`\beta`, and 
    :math:`\Gamma(z) = \int_0^\infty t^{z - 1}e^{-t}\,dt`.

    We generate noise on a uniform cartesian mesh
           
           >>> from fipy.variables.variable import Variable
           >>> alpha = Variable()
           >>> beta = Variable()
           >>> from fipy.meshes.grid2D import Grid2D
           >>> noise = GammaNoiseVariable(mesh = Grid2D(nx = 100, ny = 100), shape = alpha, rate = beta)
           
    We histogram the root-volume-weighted noise distribution
    
           >>> from fipy.variables.histogramVariable import HistogramVariable
           >>> histogram = HistogramVariable(distribution = noise, dx = 0.1, nx = 300)
           
    and compare to a Gaussian distribution
    
           >>> from fipy.variables.cellVariable import CellVariable
           >>> gammadist = CellVariable(mesh = histogram.getMesh())
           >>> x = histogram.getMesh().getCellCenters()[0]
           
           >>> if __name__ == '__main__':
           ...     from fipy import Viewer
           ...     viewer = Viewer(vars=noise, datamin=0, datamax=30)
           ...     histoplot = Viewer(vars=(histogram, gammadist), 
           ...                        datamin=0, datamax=1)
           
           >>> from fipy.tools.numerix import arange, exp
           >>> from scipy.special import gamma as Gamma
           
           >>> for shape in arange(1,8,1):
           ...     alpha.setValue(shape)
           ...     for rate in arange(0.5,2.5,0.5):
           ...         beta.setValue(rate)
           ...         gammadist.setValue(x**(alpha - 1) * (beta**alpha * exp(-beta * x)) / Gamma(alpha))
           ...         if __name__ == '__main__':
           ...             import sys
           ...             print >>sys.stderr, "alpha: %g, beta: %g" % (alpha, beta)
           ...             viewer.plot()
           ...             histoplot.plot()

           >>> print abs(noise.getFaceGrad().getDivergence().getCellVolumeAverage()) < 5e-15
           1

    .. image:: fipy/variables/gamma.*
      :scale: 25
      :align: center

    .. image:: fipy/variables/gamma-histogram.*
      :scale: 25
      :align: center

    """
    def __init__(self, mesh, shape, rate, name = '', hasOld = 0):
        r"""
        :Parameters:
            - `mesh`: The mesh on which to define the noise.
            - `shape`: The shape parameter, :math:`\alpha`.
            - `rate`: The rate or inverse scale parameter, :math:`\beta`.
                 
        """
        NoiseVariable.__init__(self, mesh = mesh, name = name, hasOld = hasOld)
        self.shapeParam = self._requires(shape)
        self.rate = self._requires(rate)
    
    def random(self):
        return random.gamma(shape=self.shapeParam, scale=self.rate, 
                            size=[self.getMesh().globalNumberOfCells])

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 

