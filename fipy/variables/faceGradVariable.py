#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "faceGradVariable.py"
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.variables.faceVariable import FaceVariable
from fipy.tools import numerix
from fipy.tools import inline

class _FaceGradVariable(FaceVariable):
    """
    Test case for a vector cell variable

    >>> from fipy import *
    >>> m = Grid2D(nx=3, ny=3)
    >>> x, y = m.cellCenters
    >>> v = CellVariable(mesh=m, elementshape=(3,))
    >>> v[0] = x
    >>> v[1] = y
    >>> v[2] = x**2
    >>> v0 = CellVariable(mesh=m, value=x)
    >>> v1 = CellVariable(mesh=m, value=y)
    >>> v2 = CellVariable(mesh=m, value=x**2)
    >>> v.faceGrad.globalValue.shape
    (2, 3, 24)
    >>> print v0.faceGrad # doctest: +NORMALIZE_WHITESPACE
    [[ 0.5  1.   0.5  0.5  1.   0.5  0.5  1.   0.5  0.5  1.   0.5  0.   1.   1.
       0.   0.   1.   1.   0.   0.   1.   1.   0. ]
     [ 0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.   0.
       0.   0.   0.   0.   0.   0.   0.   0.   0. ]]
    >>> print (v0.faceGrad.globalValue == v.faceGrad.globalValue[:,0]).all()
    True
    >>> print (v1.faceGrad.globalValue  == v.faceGrad.globalValue[:,1]).all()
    True
    >>> print (v2.faceGrad.globalValue  == v.faceGrad.globalValue[:,2]).all()
    True

    """
    def __init__(self, var):
        FaceVariable.__init__(self, mesh=var.mesh, elementshape=(var.mesh.dim,) + var.shape[:-1])
        self.var = self._requires(var)

    def _calcValue(self):
        if inline.doInline and self.var.rank == 0:
            return self._calcValueInline()
        else:
            return self._calcValueNoInline()

    def _calcValueInline(self):

        id1, id2 = self.mesh._adjacentCellIDs

        tangents1 = self.mesh._faceTangents1
        tangents2 = self.mesh._faceTangents2

        val = self._array.copy()

        faceNormals = self.mesh._orientedFaceNormals
        if numerix.MA.isMaskedArray(faceNormals):
            faceNormals = faceNormals.filled()

        inline._runIterateElementInline("""
            int j;
            double t1grad1, t1grad2, t2grad1, t2grad2, N, N2;
            int ID1 = ITEM(id1, i, NULL);
            int ID2 = ITEM(id2, i, NULL);

            if ITEM(exteriorFaces, i, NULL) {
                 N2 = ITEM(facevar, i, NULL);
            } else {
                 N2 = ITEM(var, ID2, NULL);
            }

            N = (N2 - ITEM(var, ID1, NULL)) / ITEM(dAP, i, NULL);

            t1grad1 = t1grad2 = t2grad1 = t2grad2 = 0.;

            t1grad1 += ITEM(tangents1, i, vec) * ITEM(cellGrad, ID1, vec);
            t1grad2 += ITEM(tangents1, i, vec) * ITEM(cellGrad, ID2, vec);
            t2grad1 += ITEM(tangents2, i, vec) * ITEM(cellGrad, ID1, vec);
            t2grad2 += ITEM(tangents2, i, vec) * ITEM(cellGrad, ID2, vec);

            ITEM(val, i, vec) =  ITEM(normals, i, vec) * N;
            ITEM(val, i, vec) += ITEM(tangents1, i, vec) * (t1grad1 + t1grad2) / 2.;
            ITEM(val, i, vec) += ITEM(tangents2, i, vec) * (t2grad1 + t2grad2) / 2.;
        """,tangents1 = tangents1,
            tangents2 = tangents2,
            cellGrad = self.var.grad.numericValue,
            normals = faceNormals,
            id1 = id1,
            id2 = id2,
            dAP = numerix.array(self.mesh._cellDistances),
            var = self.var.numericValue,
            facevar = self.var.faceValue.numericValue,
            exteriorFaces = self.mesh.exteriorFaces.numericValue,
            val = val,
            ni = tangents1.shape[1],
            shape=numerix.array(numerix.shape(tangents1)))

        return self._makeValue(value = val)

    def _calcValueNoInline(self):
        dAP = self.mesh._cellDistances
        id1, id2 = self.mesh._adjacentCellIDs

        N2 = numerix.take(self.var.value,id2, axis=-1)

        faceMask = numerix.array(self.mesh.exteriorFaces)

        ## The following conditional is required because empty
        ## indexing is not altogether functional.  This
        ## numpy.empty((0,))[[]] and this numpy.empty((0,))[...,[]]
        ## both work, but this numpy.empty((3, 0))[...,[]] is
        ## broken.

        if self.var.faceValue.shape[-1] != 0:
            s = (Ellipsis, faceMask)
        else:
            s = (faceMask,)

        N2[s] = self.var.faceValue[s]

        N = (N2 - numerix.take(self.var,id1, axis=-1)) / dAP

        normals = self.mesh._orientedFaceNormals

        tangents1 = self.mesh._faceTangents1
        tangents2 = self.mesh._faceTangents2
        cellGrad = self.var.grad.numericValue

        grad1 = numerix.take(cellGrad, id1, axis=-1)
        grad2 = numerix.take(cellGrad, id2, axis=-1)

        s = (slice(0,None,None),) + (numerix.newaxis,) * (len(grad1.shape) - 2) + (slice(0,None,None),)
        t1grad1 = numerix.sum(tangents1[s] * grad1, 0)
        t1grad2 = numerix.sum(tangents1[s] * grad2, 0)
        t2grad1 = numerix.sum(tangents2[s] * grad1, 0)
        t2grad2 = numerix.sum(tangents2[s] * grad2, 0)

        T1 = (t1grad1 + t1grad2) / 2.
        T2 = (t2grad1 + t2grad2) / 2.

        return normals[s] * N[numerix.newaxis] + tangents1[s] * T1[numerix.newaxis] + tangents2[s] * T2[numerix.newaxis]

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
