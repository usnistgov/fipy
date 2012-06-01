#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "abstractTopology.py"
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
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix

__all__ = []

class _AbstractTopology(object):
    _concatenatedClass = None

    def __init__(self, mesh):
        self.mesh = mesh
        
    @property
    def _isOrthogonal(self):
        raise NotImplementedError

    @property
    def _globalNonOverlappingCellIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh. 
        
        Does not include the IDs of boundary cells.

        E.g., would return [0, 1, 4, 5] for mesh A

            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfCells)

    @property
    def _globalOverlappingCellIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh. 
        
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 4, 5, 6] for mesh A

            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfCells)

    @property
    def _localNonOverlappingCellIDs(self):
        """Return the IDs of the local mesh in isolation. 
        
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

            A        B
        ------------------
        | 3 | 4 || 4 | 5 |
        ------------------
        | 0 | 1 || 1 | 2 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfCells)

    @property
    def _localOverlappingCellIDs(self):
        """Return the IDs of the local mesh in isolation. 
        
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3, 4, 5] for mesh A

            A        B
        ------------------
        | 3 | 4 || 5 |   |
        ------------------
        | 0 | 1 || 2 |   |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfCells)

    @property
    def _globalNonOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh. 
        
        Does not include the IDs of faces of boundary cells.

        E.g., would return [0, 1, 4, 5, 8, 9, 12, 13, 14, 17, 18, 19]
        for mesh A

            A   ||   B
        --8---9---10--11--
       17   18  19  20   21
        --4---5----6---7--
       12   13  14  15   16
        --0---1----2---3--
                ||
                
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfFaces)

    @property
    def _globalOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh. 
        
        Includes the IDs of faces of boundary cells.
        
        E.g., would return [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 
        14, 15, 17, 18, 19, 20] for mesh A

            A   ||   B
        --8---9---10--11--
       17   18  19  20   21
        --4---5----6---7--
       12   13  14  15   16
        --0---1----2---3--
                ||
                
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfFaces)

    @property
    def _localNonOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in isolation. 
        
        Does not include the IDs of faces of boundary cells.
        
        E.g., would return [0, 1, 3, 4, 6, 7, 9, 10, 11, 13, 14, 15]
        for mesh A

            A   ||   B
        --6---7-----7---8--
       13   14 15/14 15   16
        --3---4-----4---5--
        9   10 11/10 11   12
        --0---1-----1---2--
                ||
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfFaces)

    @property
    def _localOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in isolation. 
        
        Includes the IDs of faces of boundary cells.
        
        E.g., would return [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
        12, 13, 14, 15, 16] for mesh A

            A   ||   B
        --6---7----8------
       13   14  15  16   |
        --3---4----5------
        9   10  11  12   |
        --0---1----2------
                ||
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.numberOfFaces)
     
    # abstract element types mutually understood by FiPy and other meshing systems
    # (Vtk, Gmsh, etc.)
    _elementTopology = dict([(k, v) for (v, k) in enumerate(("vertex",
                                                             "line",
                                                             "triangle", 
                                                             "quadrangle", 
                                                             "pixel",
                                                             "polygon",
                                                             "tetrahedron", 
                                                             "hexahedron", 
                                                             "voxel",
                                                             "prism", 
                                                             "pyramid", 
                                                             "unknown"))])

    @property
    def _cellTopology(self):
        """return a map of the topology of each cell of grid"""
        raise NotImplementedError

