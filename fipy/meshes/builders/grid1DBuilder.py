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

from fipy.meshes.builders.abstractGridBuilder import _AbstractGridBuilder

from fipy.meshes.builders.utilityClasses import (_UniformNumPts,
                                                 _NonuniformNumPts,
                                                 _DOffsets,
                                                 _UniformOrigin)
from fipy.tools import numerix
 
class _Grid1DBuilder(_AbstractGridBuilder):

    def buildGridData(self, *args, **kwargs):
        kwargs["cacheOccupiedNodes"] = True
        super(_Grid1DBuilder, self).buildGridData(*args, **kwargs)

    def _packOverlap(self, first, second):
        return {'left': first, 'right': second}

    def _packOffset(self, arg):
        return arg

    @property
    def _specificGridData(self):
        return [self.occupiedNodes]

    def _calcShape(self):
        return (self.ns[0],)
             
    def _calcPhysicalShape(self):
        """Return physical dimensions of Grid1D."""
        from fipy.tools.dimensions.physicalField import PhysicalField
        return PhysicalField(value = (self.ns[0] * self.ds[0] * self.scale,))
                      
    def _calcMeshSpacing(self):
        return numerix.array((self.ds[0],))[...,numerix.newaxis]

    @staticmethod
    def createVertices(dx, nx):
        x = _AbstractGridBuilder.calcVertexCoordinates(dx, nx)
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
      
class _NonuniformGrid1DBuilder(_Grid1DBuilder):

    def __init__(self):
        self.NumPtsCalcClass = _NonuniformNumPts

        super(_NonuniformGrid1DBuilder, self).__init__()
 
    def buildGridData(self, *args, **kwargs):
        # call super for side-effects
        super(_NonuniformGrid1DBuilder, self).buildGridData(*args, **kwargs)

        (self.offsets, 
         self.ds) = _DOffsets.calcDOffsets(self.ds, self.ns, self.offset)

        self.vertices = _Grid1DBuilder.createVertices(self.ds[0], self.ns[0]) \
                         + ((self.offsets[0],),) 
        self.faces = _Grid1DBuilder.createFaces(self.numberOfVertices)
        self.numberOfFaces = len(self.faces[0])
        self.cells = _Grid1DBuilder.createCells(self.ns[0])

    @property
    def _specificGridData(self):
        return super(_NonuniformGrid1DBuilder, self)._specificGridData \
                + [self.vertices,
                   self.faces,
                   self.cells]
                                     
class _UniformGrid1DBuilder(_Grid1DBuilder):

    def __init__(self):
        self.NumPtsCalcClass = _UniformNumPts

        super(_UniformGrid1DBuilder, self).__init__()

    def buildGridData(self, ns, ds, overlap, communicator, origin):
        super(_UniformGrid1DBuilder, self).buildGridData(ns, ds, overlap,
                                                        communicator)

        self.origin = _UniformOrigin.calcOrigin(origin, 
                                                self.offset, self.ds, self.scale)

        if 0 in self.ns:
            self.numberOfFaces = 0
        else:
            self.numberOfFaces = self.ns[0] + 1

        self.numberOfCells = self.ns[0]

    @property
    def _specificGridData(self):
        return super(_UniformGrid1DBuilder, self)._specificGridData \
                + [self.origin]




