#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "SkewedGrid2D.py"
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

__docformat__ = 'restructuredtext'
 
from fipy.tools import numerix
from fipy.tools.numerix import random
from fipy.tools.dimensions.physicalField import PhysicalField

from fipy.meshes.mesh2D import Mesh2D
from fipy.meshes import Grid2D

__all__ = ["SkewedGrid2D"]

class SkewedGrid2D(Mesh2D):
    """
    Creates a 2D grid mesh with horizontal faces numbered first and then
    vertical faces.  The points are skewed by a random amount (between `rand`
    and `-rand`) in the X and Y directions.
    """
    def __init__(self, dx = 1., dy = 1., nx = None, ny = 1, rand = 0, *args, **kwargs):
        self.args = {
            'dx': dx, 
            'dy': dy, 
            'nx': nx, 
            'ny': ny, 
            'rand': rand
        }
        
        self.nx = nx
        self.ny = ny
        
        self.dx = PhysicalField(value = dx)
        scale = PhysicalField(value = 1, unit = self.dx.unit)
        self.dx /= scale
        
        self.dy = PhysicalField(value = dy)
        if self.dy.unit.isDimensionless():
            self.dy = dy
        else:
            self.dy /= scale

        self.grid = Grid2D(nx=nx, ny=ny, dx=dx, dy=dy)

        self.numberOfVertices = self.grid._numberOfVertices
        
        vertices = self.grid.vertexCoords

        changedVertices = numerix.zeros(vertices.shape, 'd')

        for i in range(len(vertices[0])):
            if((i % (nx+1)) != 0 and (i % (nx+1)) != nx and (i / nx+1) != 0 and (i / nx+1) != ny):
                changedVertices[0, i] = vertices[0, i] + (rand * ((random.random() * 2) - 1))
                changedVertices[1, i] = vertices[1, i] + (rand * ((random.random() * 2) - 1))
            else:
                changedVertices[0, i] = vertices[0, i]
                changedVertices[1, i] = vertices[1, i]


        faces = self.grid.faceVertexIDs
        
        cells = self.grid.cellFaceIDs

        Mesh2D.__init__(self, changedVertices, faces, cells, *args, **kwargs)
        
        self.scale = scale
        
    @property
    def physicalShape(self):
        """Return physical dimensions of Grid2D.
        """
        return PhysicalField(value = (self.nx * self.dx * self.scale, self.ny * self.dy * self.scale))

    @property
    def _meshSpacing(self):
        return numerix.array((self.dx,self.dy))[...,numerix.newaxis]
    
    @property
    def shape(self):
        return (self.nx, self.ny)
    
