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

from abstractGridBuilder import AbstractGridBuilder

from fipy.meshes.builders.utilityClasses import (UniformNumPts,
                                                 NonuniformNumPts,
                                                 DOffsets,
                                                 UniformOrigin)
from fipy.tools import numerix
 
class Grid1DBuilder(AbstractGridBuilder):

    def getParallelInfo(self):
        """
        :Returns:
            - A tuple.
        """
        return list(super(Grid1DBuilder, self).getParallelInfo()) \
                + [self.occupiedNodes]

    def _packOverlap(self, first, second):
        return {'left': first, 'right': second}

    def _packOffset(self, arg):
        return arg
    
    @staticmethod
    def createVertices(dx, nx):
        x = AbstractGridBuilder.calcVertexCoordinates(dx, nx)
        return x[numerix.newaxis,...]
    
    @staticmethod
    def createFaces(numVertices):
        if numVertices == 1:
            return numerix.arange(0)[numerix.newaxis, ...]
        else:
            return numerix.arange(numVertices)[numerix.newaxis, ...]

    @staticmethod
    def createCells(nx):
        """
        cells = (f1, f2) going left to right.
        f1 etc. refer to the faces
        """
        f1 = numerix.arange(nx)
        f2 = f1 + 1
        return numerix.array((f1, f2))
      
class NonuniformGrid1DBuilder(Grid1DBuilder):

    def __init__(self):
        self.NumPtsCalcClass = NonuniformNumPts

        super(NonuniformGrid1DBuilder, self).__init__()
 
    def buildPostParallelGridInfo(self, ns, ds, offset):
        # call super for side-effects
        super(NonuniformGrid1DBuilder, self).buildPostParallelGridInfo(ns)

        (self.offsets, 
         self.newDs) = DOffsets.calcDOffsets(ds, ns, offset)

        self.vertices = Grid1DBuilder.createVertices(ds[0], ns[0]) \
                         + ((self.offsets[0],),) 
        self.faces = Grid1DBuilder.createFaces(self.numberOfVertices)
        self.numberOfFaces = len(self.faces[0])
        self.cells = Grid1DBuilder.createCells(ns[0])

    def getPostParallelGridInfo(self):
        return (self.ns,
                self.newDs,
                self.offsets,
                self.vertices,
                self.faces,
                self.cells,
                self.numberOfVertices,
                self.numberOfFaces,
                self.numberOfCells)
                                     
class UniformGrid1DBuilder(Grid1DBuilder):

    def __init__(self):
        self.NumPtsCalcClass = UniformNumPts

        super(UniformGrid1DBuilder, self).__init__()

    def buildPostParallelGridInfo(self, ns, ds, offset, origin, scale):
        super(UniformGrid1DBuilder, self).buildPostParallelGridInfo(ns)

        self.origin = UniformOrigin.calcOrigin(origin, offset, ds, scale)

        if 0 in ns:
            self.numberOfFaces = 0
        else:
            self.numberOfFaces = ns[0] + 1

        self.numberOfCells = ns[0]

    def getPostParallelGridInfo(self):

        return (self.origin,
                self.numberOfVertices,
                self.numberOfFaces,
                self.numberOfCells)




