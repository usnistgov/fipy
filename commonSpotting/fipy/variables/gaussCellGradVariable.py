#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gaussCellGradVariable.py"
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
 
from fipy.variables.cellVariable import CellVariable
from fipy.tools import numerix
from fipy.tools import inline
from fipy.variables.faceGradContributionsVariable import _FaceGradContributions


class _GaussCellGradVariable(CellVariable):
    def __init__(self, var, name=''):
        CellVariable.__init__(self, mesh=var.getMesh(), name=name, rank=var.getRank() + 1)
        self.var = self._requires(var)
        self.faceGradientContributions = _FaceGradContributions(self.var)
        
    def _calcValueIn(self, N, M, ids, orientations, volumes):
        val = self._getArray().copy()

        inline._runIterateElementInline("""
            ITEM(val, i, vec) = 0.;

            int k;
            for (k = 0; k < M; k++) {
                int id = ITEM(ids, i, &k);
                ITEM(val, i, vec) += ITEM(orientations, i, &k) * ITEM(areaProj, id, vec) * ITEM(faceValues, id, NULL);
            }
                
            ITEM(val, i, vec) /= ITEM(volumes, i, NULL);
        """,val = val,
            ids = numerix.array(numerix.MA.filled(ids, 0)),
            orientations = numerix.array(numerix.MA.filled(orientations, 0)),
            volumes = numerix.array(volumes),
            areaProj = numerix.array(self.mesh._getAreaProjections()),
            faceValues = numerix.array(self.var.getArithmeticFaceValue()),
            M = M,
            ni = N, 
            shape=numerix.array(numerix.shape(val)))

        return self._makeValue(value = val)
            
    def _calcValuePy(self, N, M, ids, orientations, volumes):
        contributions = numerix.take(self.faceGradientContributions, ids, axis=1)

        grad = numerix.array(numerix.sum(orientations * contributions, 1))

        return grad / volumes

    def _calcValue(self):
        N = self.mesh.getNumberOfCells()
        M = self.mesh._getMaxFacesPerCell()
        
        ids = self.mesh._getCellFaceIDs()

        orientations = self.mesh._getCellFaceOrientations()
        volumes = self.mesh.getCellVolumes()

        return inline._optionalInline(self._calcValueIn, self._calcValuePy, N, M, ids, orientations, volumes)
