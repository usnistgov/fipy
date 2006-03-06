#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "periodicGrid1D.py"
 #                                    created: 11/10/03 {3:30:42 PM} 
 #                                last update: 3/4/06 {12:25:15 AM} 
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

"""
Peridoic 1D Mesh
"""
__docformat__ = 'restructuredtext'

from fipy.meshes.numMesh.grid1D import Grid1D

class PeriodicGrid1D(Grid1D):
    """
    
    Creates a Periodic grid mesh.
        
        >>> mesh = PeriodicGrid1D(dx = (1, 2, 3))
        
        >>> print mesh.getExteriorFaces()
        [3,]

        >>> print mesh.getFaceCellIDs()
        [[2 ,0 ,]
         [0 ,1 ,]
         [1 ,2 ,]
         [2 ,-- ,]]

        >>> print mesh._getCellDistances()
        [ 2. , 1.5, 2.5, 1.5,] 1

        >>> print mesh._getCellToCellDistances()
        [[ 2. , 1.5,]
         [ 1.5, 2.5,]
         [ 2.5, 2. ,]] 1

        >>> print mesh._getFaceNormals()
        [[ 1.,]
         [ 1.,]
         [ 1.,]
         [ 1.,]]

        >>> print mesh._getCellVertexIDs()
        [[1,0,]
         [2,1,]
         [2,0,]]
    """
    def __init__(self, dx = 1., nx = None):
        Grid1D.__init__(self, dx = dx, nx = nx)
        self._connectFaces(self.getFacesLeft(), self.getFacesRight())

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
