from fipy.boundaryConditions import FixedValue
from fipy.tools.numerix import take, where
from fipy.variables import FaceVariable

class _FixedBCFaceGradVariable(FaceVariable):
    def __init__(self, var, boundaryConditions=()):
        FaceVariable.__init__(self, mesh=var.getMesh(), rank=var.getRank() + 1)
        self.var = self._requires(var)
        self.bcs = boundaryConditions

    def _calcValue(self):        
        from fipy.tools import inline
        return inline._optionalInline(self._calcValueInline, self._calcValuePy)
    
    def _calcValuePy(self):
        dAP = self.mesh._getCellDistances().getValue()
        id1, id2 = [id.getValue() for id in self.mesh._getAdjacentCellIDs()]
        v1 = take(self.var.getValue(), id1)
        v2 = take(self.var.getValue(), id2)
        
        for bc in self.bcs:
            if isinstance(bc, FixedValue):
                v2[bc.faces.getValue()] = bc._getValue()
        
        N = (v2 - v1) / dAP
        normals = self.mesh._getOrientedFaceNormals().getValue()
        
        tangents1 = self.mesh._getFaceTangents1().getValue()
        tangents2 = self.mesh._getFaceTangents2().getValue()
        cellGrad = self.var.getGrad().getValue()
        
        grad1 = take(cellGrad, id1, axis=1)
        grad2 = take(cellGrad, id2, axis=1)
        
        t1grad1 = sum(tangents1*grad1,0)
        t1grad2 = sum(tangents1*grad2,0)
        t2grad1 = sum(tangents2*grad1,0)
        t2grad2 = sum(tangents2*grad2,0)
        
        T1 = (t1grad1 + t1grad2) / 2.
        T2 = (t2grad1 + t2grad2) / 2.

        T1 = (id1 == id2) * T1
        T2 = (id1 == id2) * T2

        return normals * N + tangents1 * T1 + tangents2 * T2

    def _calcValueInline(self):
        return self._calcValuePy()
