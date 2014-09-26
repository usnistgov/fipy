#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cylindricalNonUniformGrid1D.py"
 #
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
 #  
 # ###################################################################
 ##

"""
1D Mesh
"""
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import parallelComm

from fipy.meshes.nonUniformGrid1D import NonUniformGrid1D

__all__ = ["CylindricalNonUniformGrid1D"]

class CylindricalNonUniformGrid1D(NonUniformGrid1D):
    """
    Creates a 1D cylindrical grid mesh.
    
        >>> mesh = CylindricalNonUniformGrid1D(nx = 3)
        >>> print mesh.cellCenters
        [[ 0.5  1.5  2.5]]
         
        >>> mesh = CylindricalNonUniformGrid1D(dx = (1, 2, 3))
        >>> print mesh.cellCenters
        [[ 0.5  2.   4.5]]

        >>> print numerix.allclose(mesh.cellVolumes, (0.5, 4., 13.5)) # doctest: +PROCESSOR_0
        True
         
        >>> mesh = CylindricalNonUniformGrid1D(nx = 2, dx = (1, 2, 3))
        Traceback (most recent call last):
        ...
        IndexError: nx != len(dx)

        >>> mesh = CylindricalNonUniformGrid1D(nx=2, dx=(1., 2.)) + ((1.,),)
        >>> print mesh.cellCenters
        [[ 1.5  3. ]]
        >>> print numerix.allclose(mesh.cellVolumes, (1.5, 6)) # doctest: +PROCESSOR_0
        True
        
    """
    def __init__(self, dx=1., nx=None, origin=(0,), overlap=2, communicator=parallelComm, *args, **kwargs):
        scale = PhysicalField(value=1, unit=PhysicalField(value=dx).unit)
        self.origin = PhysicalField(value=origin)
        self.origin /= scale
    
        super(CylindricalNonUniformGrid1D, self).__init__(dx=dx,
                                                          nx=nx,
                                                          overlap=overlap,
                                                          communicator=communicator,
                                                          *args,
                                                          **kwargs)

        self.vertexCoords += origin
        self.args['origin'] = origin

    def _calcFaceCenters(self):
        faceCenters = super(CylindricalNonUniformGrid1D, self)._calcFaceCenters()
        return faceCenters + self.origin
    
    def _calcFaceAreas(self):
        return self._calcFaceCenters()[0]

    def _calcCellVolumes(self):
        return super(CylindricalNonUniformGrid1D, self)._calcCellVolumes() / 2.   

    def _translate(self, vector):
        return CylindricalNonUniformGrid1D(dx=self.args['dx'], nx=self.args['nx'], 
                                           origin=numerix.array(self.args['origin']) + vector,
                                           overlap=self.args['overlap'])
                                     
    def __mul__(self, factor):
        return CylindricalNonUniformGrid1D(dx=self.args['dx'] * factor, nx=self.args['nx'], 
                                           origin=numerix.array(self.args['origin']) * factor,
                                           overlap=self.args['overlap'])

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected. Fixed a bug where the following throws
        an error on solve() when nx is a float.

            >>> from fipy import CellVariable, DiffusionTerm
            >>> mesh = CylindricalNonUniformGrid1D(nx=3., dx=(1., 2., 3.))
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var)

        This test is for https://github.com/usnistgov/fipy/issues/372. Cell
        volumes were being returned as binOps rather than arrays.

            >>> m = CylindricalNonUniformGrid1D(dx=(1., 2., 3., 4.), nx=4)
            >>> print isinstance(m.cellVolumes, numerix.ndarray)
            True
            >>> print isinstance(m._faceAreas, numerix.ndarray)
            True

        If the above types aren't correct, the divergence operator's value can be a binOp

            >>> print isinstance(CellVariable(mesh=m).arithmeticFaceValue.divergence.value, numerix.ndarray)
            True
            
        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
