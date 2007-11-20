## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "faceToVertexVariable.py"
 #                                     created: 11/17/07 {9:54:32 AM}
 #                                 last update: 11/17/07 {11:04:21 AM}
 # Author: Jonathan Guyer
 # E-mail: <jguyer@his.com>
 #   mail: Alpha Cabal
 #    www: <http://alphatcl.sourceforge.net>
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
 # History
 # 
 # modified   by  rev reason
 # ---------- --- --- -----------
 # 2007-11-17 JEG 1.0 original
 # 
 # #############################################################################
 ##

from fipy.variables.vertexVariable import _VertexVariable

class _FaceToVertexVariable(_VertexVariable):
    def __init__(self, var):
        _VertexVariable.__init__(self, mesh=var.getMesh(), elementshape=var.shape[:-1])
        self.var = self._requires(var)

    def _calcValue(self):
        from fipy.tools import numerix, vector
        
        V = self.mesh._getNumberOfVertices()
        value = numerix.zeros(self.var.shape[:-1] + (V,), 'd')
        weight = numerix.zeros((V,), 'l')
        M, F = self.mesh._getFaceVertexIDs().shape
        mask = numerix.zeros((F,), 'l')

        for j in range(M):
            ids = self.mesh._getFaceVertexIDs()[j]
            mask[:] = ids.getMask()
            vector.putAdd(value, ids, self.var)
            vector.putAdd(weight, ids, mask == False)
                
        return value / weight

        
