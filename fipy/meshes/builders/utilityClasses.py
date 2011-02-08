#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
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
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import numerix

class UniformOrigin(object):
    """
    Used to calculate the origin for uniform grids in a
    dimensionally-independent way.
    """

    @staticmethod
    def calcOrigin(origin, offset, ds, scale):
        newOrigin  = PhysicalField(value=origin)
        newOrigin /= scale

        if type(offset) in [int, float]:
            newOrigin += offset * ds[0]
        else:
            newOrigin += (o*d for o, d in zip(offset, ds))

        return newOrigin

class DOffsets(object):
    """
    For use by non-uniform grid builders.
    """

    @staticmethod
    def calcDOffsets(ds, ns, offset):
        """
        :Parameters:
            - `ds`: A list, e.g. [dx, dy]
            - `ns`: A list, e.g. [nx, ny, nz]
            - `offset`: A scalar 

        :Returns:
            - `offsetList`: a list which contains the analogous scalars to
              `XOffset`, `YOffset`, and `ZOffset`, whichever are applicable for
              the dimensionality.
            - `newDs`: a list containing proper [dx, [dy, ...]] values
        """
        offsetList = []
        newDs = []
        
        for d, n in zip(ds, ns):
            if numerix.getShape(d) is not ():
                offsetList.append(numerix.sum(d[0:offset]))
                newDs.append(d[offset:offset + n])
            else:
                offsetList.append(d * offset)
                newDs.append(d)

        return offsetList, newDs

class AbstractNumPts(object):
    """
    Interface definition for NumPtsCalculators.
    """

    @staticmethod
    def calcNs(ns, ds):
        raise NotImplementedError
 
class NonuniformNumPts(AbstractNumPts):
    """
    For use by non-uniform grid builders.
    """

    @staticmethod
    def calcNs(ns, ds):
        axis = ["x", "y", "z"][:len(ns)]
        newNs = []

        for a, d, n in zip(axis, ds, ns):
            newNs.append(NonuniformNumPts._calcNumPts(d=d, n=n, axis=a))

        return newNs

    @staticmethod
    def _calcNumPts(d, n = None, axis = "x"):
        """
        Calculate the number of cells along the specified axis, based
        on either the specified number or on the number elements in the
        cell  `d` spacings.
        
        Used by the `Grid` meshes.

        This tests a bug that was occuring with PeriodicGrid1D when
        using a numpy float as the argument for the grid spacing.

           >>> from fipy.meshes.periodicGrid1D import PeriodicGrid1D
           >>> PeriodicGrid1D(nx=2, dx=numerix.float32(1.))
           PeriodicGrid1D(dx=1.0, nx=2)
        """

        if type(d) in [int, float] or not hasattr(d, '__len__'):
            n = int(n or 1)
        else:
            n = int(n or len(d))
            if n != len(d) and len(d) != 1:
                raise IndexError, "n%s != len(d%s)" % (axis, axis)
                
        return n

class UniformNumPts(AbstractNumPts):
    """
    For use by uniform grid builders.
    """
 
    @staticmethod
    def calcNs(ns, ds):
        return map(lambda x: int(x), ns)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
