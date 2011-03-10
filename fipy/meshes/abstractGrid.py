#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "abstractMesh.py"
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

from fipy.tools import serial
from fipy.tools import numerix
from fipy.tools.decorators import getsetDeprecated
from fipy.tools.numerix import MA

 
def AbstractGridFactory(parent):
    class AbstractGrid(parent):

        @getsetDeprecated
        def getPhysicalShape(self):
            return self.physicalShape

        @getsetDeprecated
        def _getMeshSpacing(self):
            return self._meshSpacing
       
        @getsetDeprecated
        def getShape(self):
            return self.shape

        def _isOrthogonal(self):
            return True
             
        def __getstate__(self):
            """
            Used internally to collect the necessary information to ``pickle`` the 
            `Grid2D` to persistent storage.
            """
            return self.args

        def __setstate__(self, dict):
            """
            Used internally to create a new `Grid2D` from ``pickled`` 
            persistent storage.
            """
            self.__init__(**dict)

    return AbstractGrid
                               
def AbstractGrid1DFactory(parent):

    from fipy.meshes.mesh1D import Mesh1D
    class AbstractGrid1D(AbstractGridFactory(parent)):
 
        def __repr__(self):
            if self.args["nx"] is None:
                return "%s(dx=%s)" % (self.__class__.__name__, 
                                      str(self.args["dx"]))
            else:
                return "%s(dx=%s, nx=%d)" % (self.__class__.__name__, 
                                             str(self.args["dx"]), 
                                             self.args["nx"])
 
        @property
        def _concatenatedClass(self):
            return Mesh1D
                                                                
        @property
        def _globalNonOverlappingCellIDs(self):
            """
            Return the IDs of the local mesh in the context of the
            global parallel mesh. Does not include the IDs of boundary cells.

            E.g., would return [0, 1] for mesh A

                A        B
            ------------------
            | 0 | 1 || 2 | 3 |
            ------------------
            
            .. note:: Trivial except for parallel meshes
            """

            return numerix.arange(self.offset + self.overlap['left'], 
                                  self.offset + self.nx - self.overlap['right'])

        @property
        def _globalOverlappingCellIDs(self):
            """
            Return the IDs of the local mesh in the context of the
            global parallel mesh. Includes the IDs of boundary cells.
            
            E.g., would return [0, 1, 2] for mesh A

                A        B
            ------------------
            | 0 | 1 || 2 | 3 |
            ------------------
            
            .. note:: Trivial except for parallel meshes
            """
            return numerix.arange(self.offset, self.offset + self.nx)

        @property
        def _localNonOverlappingCellIDs(self):
            """
            Return the IDs of the local mesh in isolation. 
            Does not include the IDs of boundary cells.
            
            E.g., would return [0, 1] for mesh A

                A        B
            ------------------
            | 0 | 1 || 1 | 2 |
            ------------------
            
            .. note:: Trivial except for parallel meshes
            """
            return numerix.arange(self.overlap['left'], 
                                  self.nx - self.overlap['right'])

        @property
        def _localOverlappingCellIDs(self):
            """
            Return the IDs of the local mesh in isolation. 
            Includes the IDs of boundary cells.
            
            E.g., would return [0, 1, 2] for mesh A

                A        B
            ------------------
            | 0 | 1 || 2 |   |
            ------------------
            
            .. note:: Trivial except for parallel meshes
            """
            return numerix.arange(0, self.nx)

        @property
        def _globalNonOverlappingFaceIDs(self):
            """
            Return the IDs of the local mesh in the context of the
            global parallel mesh. Does not include the IDs of boundary cells.

            E.g., would return [0, 1, 2] for mesh A

                A    ||   B
            ------------------
            0   1    2   3   4
            ------------------
            
            .. note:: Trivial except for parallel meshes
            """
            return numerix.arange(self.offset + self.overlap['left'], 
                                  self.offset + self.numberOfFaces - self.overlap['right'])

        @property
        def _globalOverlappingFaceIDs(self):
            """
            Return the IDs of the local mesh in the context of the
            global parallel mesh. Includes the IDs of boundary cells.
            
            E.g., would return [0, 1, 2, 3] for mesh A

                A    ||   B
            ------------------
            0   1    2   3   4
            ------------------
            
            .. note:: Trivial except for parallel meshes
            """
            return numerix.arange(self.offset, self.offset + self.numberOfFaces)

        @property
        def _localNonOverlappingFaceIDs(self):
            """
            Return the IDs of the local mesh in isolation. 
            Does not include the IDs of boundary cells.
            
            E.g., would return [0, 1, 2] for mesh A

                A    ||   B
            ------------------
            0   1   2/1  2   3
            ------------------
            
            .. note:: Trivial except for parallel meshes
            """
            return numerix.arange(self.overlap['left'], 
                                  self.numberOfFaces - self.overlap['right'])

        @property
        def _localOverlappingFaceIDs(self):
            """
            Return the IDs of the local mesh in isolation. 
            Includes the IDs of boundary cells.
            
            E.g., would return [0, 1, 2, 3] for mesh A

                A   ||   B
            ------------------
            0   1   2   3    |
            ------------------
            
            .. note:: Trivial except for parallel meshes
            """
            return numerix.arange(0, self.numberOfFaces)
    
    return AbstractGrid1D
     
def AbstractGrid2DFactory(parent):

    from fipy.meshes.mesh2D import Mesh2D

    class AbstractGrid2D(AbstractGridFactory(parent)):
    
        def __repr__(self):
            return "%s(dx=%s, dy=%s, nx=%s, ny=%s)" \
                % (self.__class__.__name__, str(self.args["dx"]), str(self.args["dy"]), 
                   str(self.args["nx"]), str(self.args["ny"]))
                 
        @property
        def _concatenatedClass(self):
            return Mesh2D
         
        @property
        def _globalNonOverlappingCellIDs(self):
            """
            Return the IDs of the local mesh in the context of the
            global parallel mesh. Does not include the IDs of boundary cells.

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
            return numerix.arange((self.offset[1] + self.overlap['bottom']) * self.nx, 
                                  (self.offset[1] + self.ny - self.overlap['top']) * self.nx)

        @property
        def _globalOverlappingCellIDs(self):
            """
            Return the IDs of the local mesh in the context of the
            global parallel mesh. Includes the IDs of boundary cells.
            
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
            return numerix.arange(self.offset[1] * self.nx, (self.offset[1] + self.ny) * self.nx)

        @property
        def _localNonOverlappingCellIDs(self):
            """
            Return the IDs of the local mesh in isolation. 
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
            return numerix.arange(self.overlap['bottom'] * self.nx, 
                                  (self.ny - self.overlap['top']) * self.nx)

        @property
        def _localOverlappingCellIDs(self):
            """
            Return the IDs of the local mesh in isolation. 
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
            return numerix.arange(0, self.ny * self.nx)
         
    return AbstractGrid2D

def AbstractGrid3DFactory(parent):

    from fipy.meshes.mesh import Mesh

    class AbstractGrid3D(AbstractGridFactory(parent)):
     
        def __repr__(self):
            return "%s(dx=%s, dy=%s, dz=%s, nx=%d, ny=%d, nz=%d)" \
                % (self.__class__.__name__, str(self.args["dx"]), str(self.args["dy"]), str(self.args["dz"]), 
                   self.args["nx"], self.args["ny"], self.args["nz"])
     
        @property
        def _concatenatedClass(self):
            return Mesh
         
        @property
        def _globalNonOverlappingCellIDs(self):
            """
            Return the IDs of the local mesh in the context of the
            global parallel mesh. Does not include the IDs of boundary cells.
            
            .. note:: Trivial except for parallel meshes
            """
            return numerix.arange((self.offset[2] + self.overlap['front']) * self.nx * self.ny, 
                                  (self.offset[2] + self.nz - self.overlap['back']) * self.nx * self.ny)

        @property
        def _globalOverlappingCellIDs(self):
            """
            Return the IDs of the local mesh in the context of the
            global parallel mesh. Includes the IDs of boundary cells.
            
            .. note:: Trivial except for parallel meshes
            """
            
            return numerix.arange(self.offset[2] * self.nx * self.ny, (self.offset[2] + self.nz) * self.nx * self.ny)

        @property
        def _localNonOverlappingCellIDs(self):
            """
            Return the IDs of the local mesh in isolation. 
            Does not include the IDs of boundary cells.
            
            .. note:: Trivial except for parallel meshes
            """
            return numerix.arange(self.overlap['front'] * self.nx * self.ny, 
                                  (self.nz - self.overlap['back']) * self.nx * self.ny)

        @property
        def _localOverlappingCellIDs(self):
            """
            Return the IDs of the local mesh in isolation. 
            Includes the IDs of boundary cells.
            
            .. note:: Trivial except for parallel meshes
            """
            return numerix.arange(0, self.ny * self.nx * self.nz)
         
    return AbstractGrid3D
