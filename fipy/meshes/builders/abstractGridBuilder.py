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

from itertools import permutations

from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import numerix
 
class AbstractGridBuilder(object):

    def buildPreParallelGridInfo(self, ds, ns):
        """
        Build and save any information derivable from (dx, dy, dz) and 
        (nx, ny, nz). Generalized to handle any dimension. Has side-effects.

        :Parameters:
            - `ds` - A list containing grid spacing information, e.g. [dx, dy]
            - `ns` - A list containing number of grid points, e.g. [nx, ny, nz]
        """

        dim = len(ns)

        newDs = []

        newdx = PhysicalField(value = ds[0])
        scale = PhysicalField(value = 1, unit = newdx.unit)
        newdx /= scale

        newDs.append(newdx)

        for d in ds[1:]: # for remaining ds
            newD = PhysicalField(value = d)

            if newD.unit.isDimensionless():
                newD = d
            else:
                newD /= scale

            newDs.append(newD)

        newNs = self._calcNs(ns, newDs)

        globalNumCells = reduce(self._mult, newNs)
        globalNumFaces = self._calcNumFaces(newNs)

        """
        Side-effects
        """
        self.dim     = dim
        self.selfds  = newDs
        self.ns      = newNs
        self.scale   = scale

        self.globalNumberOfCells = globalNumCells
        self.globalNumberOfFaces = globalNumFaces

    def getPreParallelGridInfo(self):
        return (self.ns,
                self.selfds,
                self.dim,
                self.scale,
                self.globalNumberOfCells,
                self.globalNumberOfFaces)

    def _calcNumFaces(self, ns):
        """
        Dimensionally independent face-number calculation.

        >>> from fipy.meshes.builders import *

        >>> gb = Grid1DBuilder()
        >>> gb._calcNumFaces([1])
        2

        >>> gb2 = Grid2DBuilder()
        >>> gb2._calcNumFaces([2, 3])
        17

        >>> gb3 = Grid3DBuilder()
        >>> gb3._calcNumFaces([2, 3, 2])
        52
        """
        assert type(ns) is list

        nIter = list(permutations(range(len(ns))))

        # ensure len(nIter) == len(ns) && nIter[i][0] unique
        if len(ns) == 3: 
            nIter = nIter[::2]

        numFaces = 0
        
        for idx in nIter: 
            temp = ns[idx[0]] + 1 # base of (n_1 + 1)

            for otherIdx in idx[1:]: # for any extra dimensions beyond x
                temp *= ns[otherIdx] # build (n_x * ... + (n_1 + 1))
                                     # for each valid ordering
            numFaces += temp

        return numFaces
        
    def _calcNs(self, ns, ds):
        return self.NumPtsCalcClass.calcNs(ns, ds)

    def buildParallelInfo(self, ns, overlap, communicator):
        """
        Dimension-independent method for establishing parallel grid
        information. Has side-effects.

        Dimension specific functionality is built into `_buildOverlap` and
        `_packOffset`, which are overridden by children of this class.

        :Parameters:
            - `ns`: A tuple containing `nx`, `ny`, and `nz`, if available for
              the particular dimension of this builder. Could be `(nx)`, `(nx,
              ny)`, etc.
        """
        
        ns = list(ns)

        procID = communicator.procID
        Nproc = communicator.Nproc

        overlap = min(overlap, ns[-1])
        cellsPerNode = max(int(ns[-1] / Nproc), overlap)

        occupiedNodes = min(int(ns[-1] / (cellsPerNode or 1)), Nproc) 

        (firstOverlap,
         secOverlap,
         overlap) = self._buildOverlap(overlap, procID, occupiedNodes)

        offsetArg = min(procID, occupiedNodes-1) * cellsPerNode - firstOverlap

        local_n = cellsPerNode * (procID < occupiedNodes)

        if procID == occupiedNodes - 1:
            local_n += (ns[-1] - cellsPerNode * occupiedNodes)

        local_n += firstOverlap + secOverlap

        """
        Side-effects
        """
        self.local_ns = tuple(ns[:-1] + [local_n])
        self.occupiedNodes = occupiedNodes
        self.offset = self._packOffset(offsetArg)
        self.overlap = overlap

    def getParallelInfo(self):
        return (self.local_ns, 
                self.overlap, 
                self.offset)

    def buildPostParallelGridInfo(self, ns):
        """
        Builds final information about the grid after `buildParallelInfo` has
        been called. Establishes the follow information:

            * [nx, [ny, [nz]]]
            * numberOfHorizontalRows
            * numberOfVerticalColumns
            * numberOfLayersDeep
            * numberOfVertices
            * numberOfCells

        All other attributes needed by grids are established by children who
        override `buildPostParallelGridInfo` and call up to this incarnation.

        :Note: 
            - `spatialNums` is a list whose elements are analogous to
                * numberOfHorizontalRows
                * numberOfVerticalColumns
                * numberOfLayersDeep
              though `spatialNums` may be of length 1, 2, or 3 depending on
              dimensionality.
 
        :Parameters:
            - `ns`: A tuple containing `nx`, `ny`, and `nz`, if available for
              the particular dimension of this builder. Could be `(nx)`, `(nx,
              ny)`, etc. 
        """
        
        if 0 in ns:
            newNs = [0 for n in ns]
            spatialNums = [0 for n in ns]
        else:
            newNs = ns
            spatialNums = [n + 1 for n in ns]

        numVertices = reduce(self._mult, spatialNums)
        numCells = reduce(self._mult, newNs)

        """
        Side-effects
        """
        self.ns = newNs
        self.spatialNums = spatialNums
        self.numberOfVertices = numVertices
        self.numberOfCells = numCells

    def _buildOverlap(self, overlap, procID, occupiedNodes):
         
        (first, sec) = self._calcFirstAndSecOverlap(overlap, procID, occupiedNodes)
        return first, sec, self._packOverlap(first, sec)

    def _calcFirstAndSecOverlap(self, overlap, procID, occupiedNodes):
        """
        Only consolidated to prevent duplication in `PeriodicGrid1DBuilder`.
        """
        return (overlap * (procID > 0) * (procID < occupiedNodes),
                overlap * (procID < occupiedNodes - 1))
         
    def _packOverlap(self, first, sec):
        raise NotImplementedError

    def _packOffset(self, arg):
        raise NotImplementedError
    
    @staticmethod
    def calcVertexCoordinates(d, n):
        """
        Calculate the positions of the vertices along an axis, based on the 
        specified `Cell` `d` spacing or list of `d` spacings.
        
        Used by the `Grid` meshes.
        """
        x = numerix.zeros((n + 1), 'd')
        if n > 0:
            x[1:] = d
        return numerix.add.accumulate(x)

    def _mult(self, x, y):
        return x*y

if __name__ == '__main__':
    import doctest
    doctest.testmod()
