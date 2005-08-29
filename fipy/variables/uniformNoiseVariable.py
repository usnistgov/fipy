## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "uniformNoiseVariable.py"
 #                                    created: 8/26/05 {3:08:48 PM} 
 #                                last update: 8/26/05 {5:20:10 PM} 
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
 # protection and is in the public domain.  uniformNoiseVariable.py
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

from RandomArray import uniform

from fipy.variables.noiseVariable import NoiseVariable

class UniformNoiseVariable(NoiseVariable):
    r"""
    Represents a uniform distribution of random numbers.
    
    We generate noise on a uniform cartesian mesh
    
        >>> from fipy.meshes.grid2D import Grid2D
        >>> noise = UniformNoiseVariable(mesh = Grid2D(nx = 100, ny = 100))
        
    and histogram the noise
        
        >>> from fipy.variables.histogramVariable import HistogramVariable
        >>> histogram = HistogramVariable(distribution = noise, dx = 0.01, nx = 120, offset = -.1)
        
        >>> if __name__ == 'main':
        ...     from fipy import viewers
        ...     viewer = viewers.make(vars = noise, limits = {'datamin':0, 'datamax':1})
        ...     histoplot = viewers.make(vars = histogram)
        
        >>> for i in range(10):
        ...     noise.scramble()
        ...     if __name__ == '__main__':
        ...         viewer.plot()
        ...         histoplot.plot()

    .. image:: fipy/variables/uniform.jpg
       :scale: 25
       :align: center

    .. image:: fipy/variables/uni-histogram.pdf
       :scale: 25
       :align: center
    """
    def __init__(self, mesh, name = '', minimum = 0., maximum = 1., hasOld = 0):
        """
        :Parameters:
            - `mesh`: The mesh on which to define the noise.
            - `minimum`: The minimum (not-inclusive) value of the distribution.
            - `maximum`: The maximum (not-inclusive) value of the distribution.
        """
        self.minimum = minimum
        self.maximum = maximum
        NoiseVariable.__init__(self, mesh = mesh, name = name, hasOld = hasOld)
    
    def _calcValue(self):
        self.value = uniform(self.minimum, self.maximum,
                             shape = [self.getMesh().getNumberOfCells()])

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 

