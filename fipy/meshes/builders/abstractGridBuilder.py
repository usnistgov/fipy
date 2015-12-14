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

__all__ = []

import itertools

from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import numerix

class _AbstractGridBuilder(object):

    NumPtsCalcClass = None

    def buildGridData(self, ds, ns, overlap, communicator,
                            cacheOccupiedNodes=False):
        """
        Build and save any information relevant to the construction of a grid.
        Generalized to handle any dimension. Has side-effects.

        Dimension specific functionality is built into `_buildOverlap` and
        `_packOffset`, which are overridden by children of this class. Often,
        this method is overridden (but always called) by children classes who
        must distinguish between uniform and non-uniform behavior.

        :Note:
            - `spatialNums` is a list whose elements are analogous to
                * numberOfHorizontalRows
                * numberOfVerticalColumns
                * numberOfLayersDeep
              though `spatialNums` may be of length 1, 2, or 3 depending on
              dimensionality.

        :Parameters:
            - `ds` - A list containing grid spacing information, e.g. [dx, dy]
            - `ns` - A list containing number of grid points, e.g. [nx, ny, nz]
            - `overlap`
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
                if type(d) in [list, tuple]:
                    newD = numerix.array(d)
                else:
                    newD = d
            else:
                newD /= scale

            newDs.append(newD)


        newNs = self._calcNs(ns, newDs)

        globalNumCells = reduce(self._mult, newNs)
        globalNumFaces = self._calcGlobalNumFaces(newNs)

        """
        parallel stuff
        """

        newNs = list(newNs)

        procID = communicator.procID
        Nproc = communicator.Nproc

        overlap = min(overlap, newNs[-1])
        cellsPerNode = max(newNs[-1] // Nproc, overlap)
        occupiedNodes = min(newNs[-1] // (cellsPerNode or 1), Nproc)

        (firstOverlap,
         secOverlap,
         overlap) = self._buildOverlap(overlap, procID, occupiedNodes)

        offsetArg = min(procID, occupiedNodes-1) * cellsPerNode - firstOverlap
        offset = self._packOffset(offsetArg)

        """
        local nx, [ny, [nz]] calculation
        """
        local_n = cellsPerNode * (procID < occupiedNodes)

        if procID == occupiedNodes - 1:
            local_n += (newNs[-1] - cellsPerNode * occupiedNodes)

        local_n += firstOverlap + secOverlap

        newNs = tuple(newNs[:-1] + [local_n])

        """
        post-parallel
        """

        # "what is spatialNums?" -> see docstring
        if 0 in newNs:
            newNs = [0 for n in newNs]
            spatialNums = [0 for n in newNs]
        else:
            spatialNums = [n + 1 for n in newNs]

        spatialDict = dict(zip(["numVerticalCols",
                                "numHorizontalRows",
                                "numLayersDeep"][:len(spatialNums)],
                               spatialNums))

        numVertices = reduce(self._mult, spatialNums)
        numCells = reduce(self._mult, newNs)

        """
        Side-effects
        """
        self.dim     = dim
        self.ds      = newDs
        self.ns      = newNs
        self.scale   = scale

        self.globalNumberOfCells = globalNumCells
        self.globalNumberOfFaces = globalNumFaces

        self.offset = offset
        self.overlap = overlap

        self.spatialDict = spatialDict
        self.numberOfVertices = numVertices
        self.numberOfCells = numCells

        if cacheOccupiedNodes:
            self.occupiedNodes = occupiedNodes

    @property
    def gridData(self):
        """
        `_getSpecificGridData must be defined by children.
        """
        return self._basicGridData + self._specificGridData

    @property
    def _basicGridData(self):
        return [self.ds,
                self.ns,
                self.dim,
                self.scale,
                self.globalNumberOfCells,
                self.globalNumberOfFaces,
                self.overlap,
                self.offset,
                self.numberOfVertices,
                self.numberOfFaces,
                self.numberOfCells,
                self._calcShape(),
                self._calcPhysicalShape(),
                self._calcMeshSpacing()]

    def _calcShape(self):
        raise NotImplementedError

    def _calcPhysicalShape(self):
        raise NotImplementedError

    def _calcMeshSpacing(self):
        raise NotImplementedError

    @property
    def _specificGridData(self):
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

    def _calcGlobalNumFaces(self, ns):
        """
        Dimensionally independent face-number calculation.

        >>> from fipy.meshes.builders import *

        >>> gb = _Grid1DBuilder()
        >>> gb._calcGlobalNumFaces([1])
        2

        >>> gb2 = _Grid2DBuilder()
        >>> gb2._calcGlobalNumFaces([2, 3])
        17

        >>> gb3 = _Grid3DBuilder()
        >>> gb3._calcGlobalNumFaces([2, 3, 2])
        52
        """
        assert type(ns) is list

        # `permutations` is the cleanest way to do this, but it's new in
        # python 2.6, so we can't rely on it.
        if hasattr(itertools, "permutations"):
            nIter = list(itertools.permutations(range(len(ns))))

            # ensure len(nIter) == len(ns) && nIter[i][0] unique
            if len(ns) == 3:
                nIter = nIter[::2]
        else:
            nIter = [[[0]],
                     [[0, 1], [1, 0]],
                     [[0, 1, 2], [1, 0, 2], [2, 0, 1]]][len(ns) - 1]

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

    def _buildOverlap(self, overlap, procID, occupiedNodes):
        (first, sec) = (overlap * (procID > 0) * (procID < occupiedNodes),
                        overlap * (procID < occupiedNodes - 1))

        return first, sec, self._packOverlap(first, sec)

    def _packOverlap(self, first, sec):
        raise NotImplementedError

    def _packOffset(self, arg):
        raise NotImplementedError

    def _mult(self, x, y):
        return x*y

if __name__ == '__main__':
    import doctest
    doctest.testmod()
