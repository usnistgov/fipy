#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gridTopology.py"
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

__all__ = []

from fipy.tools import numerix

from fipy.meshes.mesh1D import Mesh1D
from fipy.meshes.mesh2D import Mesh2D
from fipy.meshes.mesh import Mesh

from fipy.meshes.topologies.abstractTopology import _AbstractTopology

class _GridTopology(_AbstractTopology):

    @property
    def _isOrthogonal(self):
        return True

class _Grid1DTopology(_GridTopology):

    _concatenatedClass = Mesh1D
                                                            
    @property
    def _globalNonOverlappingCellIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh. 
        
        Does not include the IDs of boundary cells.

        E.g., would return [0, 1] for mesh A

            A        B
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """

        return numerix.arange(self.mesh.offset + self.mesh.overlap['left'], 
                              self.mesh.offset + self.mesh.nx - self.mesh.overlap['right'])

    @property
    def _globalOverlappingCellIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh. 
        
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2] for mesh A

            A        B
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.offset, self.mesh.offset + self.mesh.nx)

    @property
    def _localNonOverlappingCellIDs(self):
        """Return the IDs of the local mesh in isolation. 
        
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1] for mesh A

            A        B
        ------------------
        | 0 | 1 || 1 | 2 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.overlap['left'], 
                              self.mesh.nx - self.mesh.overlap['right'])

    @property
    def _localOverlappingCellIDs(self):
        """Return the IDs of the local mesh in isolation. 
        
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2] for mesh A

            A        B
        ------------------
        | 0 | 1 || 2 |   |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(0, self.mesh.nx)

    @property
    def _globalNonOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh. 
        
        Does not include the IDs of faces of boundary cells.

        E.g., would return [0, 1, 2] for mesh A

            A    ||   B
        ------------------
        0   1    2   3   4
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.offset + self.mesh.overlap['left'], 
                              self.mesh.offset + self.mesh.numberOfFaces - self.mesh.overlap['right'])

    @property
    def _globalOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh. 
        
        Includes the IDs of faces of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

            A    ||   B
        ------------------
        0   1    2   3   4
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.offset, self.mesh.offset + self.mesh.numberOfFaces)

    @property
    def _localNonOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in isolation. 
        
        Does not include the IDs of faces of boundary cells.
        
        E.g., would return [0, 1, 2] for mesh A

            A    ||   B
        ------------------
        0   1   2/1  2   3
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.overlap['left'], 
                              self.mesh.numberOfFaces - self.mesh.overlap['right'])

    @property
    def _localOverlappingFaceIDs(self):
        """Return the IDs of the local mesh in isolation. 
        
        Includes the IDs of faces of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

            A   ||   B
        ------------------
        0   1   2   3    |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(0, self.mesh.numberOfFaces)
     
    @property
    def _cellTopology(self):
        """return a map of the topology of each cell of grid"""
        cellTopology = numerix.empty((self.mesh.numberOfCells,), dtype=numerix.ubyte)
        cellTopology[:] = self._elementTopology["line"]
        
        return cellTopology

class _Grid2DTopology(_GridTopology):

    _concatenatedClass = Mesh2D
     
    @property
    def _globalNonOverlappingCellIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh. 
        
        Does not include the IDs of boundary cells.

        E.g., would return [0, 1] for mesh A

        ---------
        | 4 | 5 |     
        ---------  B
        | 2 | 3 |     
        =========
        | 0 | 1 |  A
        ---------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange((self.mesh.offset[1] + self.mesh.overlap['bottom']) * self.mesh.nx, 
                              (self.mesh.offset[1] + self.mesh.ny - self.mesh.overlap['top']) * self.mesh.nx)

    @property
    def _globalOverlappingCellIDs(self):
        """Return the IDs of the local mesh in the context of the global parallel mesh. 
        
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

        ---------
        | 4 | 5 |     
        ---------  B
        | 2 | 3 |     
        =========
        | 0 | 1 |  A
        ---------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.offset[1] * self.mesh.nx, (self.mesh.offset[1] + self.mesh.ny) * self.mesh.nx)

    @property
    def _localNonOverlappingCellIDs(self):
        """Return the IDs of the local mesh in isolation. 
        
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1] for mesh A

        ---------
        | 4 | 5 |     
        ---------  B
        | 2 | 3 |     
        =========
        | 0 | 1 |  A
        ---------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.overlap['bottom'] * self.mesh.nx, 
                              (self.mesh.ny - self.mesh.overlap['top']) * self.mesh.nx)

    @property
    def _localOverlappingCellIDs(self):
        """Return the IDs of the local mesh in isolation. 
        
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

        ---------
        |   |   |     
        ---------  B
        | 2 | 3 |     
        =========
        | 0 | 1 |  A
        ---------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(0, self.mesh.ny * self.mesh.nx)
        
    @property
    def _cellTopology(self):
        """return a map of the topology of each cell of grid"""
        cellTopology = numerix.empty((self.mesh.numberOfCells,), dtype=numerix.ubyte)
        cellTopology[:] = self._elementTopology["pixel"]
        
        return cellTopology


class _Grid3DTopology(_GridTopology):
 
    _concatenatedClass = Mesh
     
    @property
    def _globalNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Does not include the IDs of boundary cells.
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange((self.mesh.offset[2] + self.mesh.overlap['front']) * self.mesh.nx * self.mesh.ny, 
                              (self.mesh.offset[2] + self.mesh.nz - self.mesh.overlap['back']) * self.mesh.nx * self.mesh.ny)

    @property
    def _globalOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Includes the IDs of boundary cells.
        
        .. note:: Trivial except for parallel meshes
        """
        
        return numerix.arange(self.mesh.offset[2] * self.mesh.nx * self.mesh.ny, 
                              (self.mesh.offset[2] + self.mesh.nz) * self.mesh.nx * self.mesh.ny)

    @property
    def _localNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Does not include the IDs of boundary cells.
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.mesh.overlap['front'] * self.mesh.nx * self.mesh.ny, 
                              (self.mesh.nz - self.mesh.overlap['back']) * self.mesh.nx * self.mesh.ny)

    @property
    def _localOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Includes the IDs of boundary cells.
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(0, self.mesh.ny * self.mesh.nx * self.mesh.nz)
     
    @property
    def _cellTopology(self):
        """return a map of the topology of each cell of grid"""
        cellTopology = numerix.empty((self.mesh.numberOfCells,), dtype=numerix.ubyte)
        cellTopology[:] = self._elementTopology["voxel"]
        
        return cellTopology

class _PeriodicGrid1DTopology(_Grid1DTopology):
    
    @property
    def _globalOverlappingCellIDs(self):
        return super(_PeriodicGrid1DTopology, self)._globalOverlappingCellIDs % self.mesh.args['nx']
