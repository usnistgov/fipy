#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "modCellGradVariable.py"
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

from fipy.variables.gaussCellGradVariable import _GaussCellGradVariable
from fipy.tools import inline
from fipy.tools import numerix

class _ModCellGradVariable(_GaussCellGradVariable):
    def __init__(self, var, modIn, modPy):
        _GaussCellGradVariable.__init__(self, var)
        self.modIn = modIn
        self.modPy = modPy


    def _calcValueInline(self, N, M, ids, orientations, volumes):
        val = self._array.copy()

        inline._runIterateElementInline(self.modIn + """
            ITEM(val, i, vec) = 0.;

            int k;
            for (k = 0; k < M; k++) {
                int id = ITEM(ids, i, &k);
                ITEM(val, i, vec) += ITEM(orientations, i, &k) * ITEM(areaProj, id, vec) * ITEM(faceValues, id, NULL);
            }

            ITEM(val, i, vec) /= ITEM(volumes, i, NULL);
            ITEM(val, i, vec) = mod(ITEM(val, i, vec) * gridSpacing[vec[0]]) /  gridSpacing[vec[0]];
        """,val = val,
            ids = numerix.array(ids),
            orientations = numerix.array(orientations),
            volumes = numerix.array(volumes),
            areaProj = numerix.array(self.mesh._areaProjections),
            faceValues = numerix.array(self.var.arithmeticFaceValue),
            M = M,
            ni = N, 
            gridSpacing = numerix.array(self.mesh._meshSpacing),
            shape=numerix.array(numerix.shape(val)))

        return self._makeValue(value = val)

    def _calcValueNoInline(self, N, M, ids, orientations, volumes):
        value = _GaussCellGradVariable._calcValueNoInline(self, N, M, ids, orientations, volumes)
        gridSpacing = self.mesh._meshSpacing
        return self.modPy(value * gridSpacing) / gridSpacing

