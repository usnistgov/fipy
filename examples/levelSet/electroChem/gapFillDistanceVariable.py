#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "gapFillDistanceVariable.py"
 #
 #  Author: Jonathan Guyer   <guyer@nist.gov>
 #  Author: Daniel Wheeler   <daniel.wheeler@nist.gov>
 #  Author: James Warren     <jwarren@nist.gov>
 #  Author: Andrew Acquaviva <andrewa@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  setup.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable
from fipy.variables.distanceVariable import DistanceVariable

class GapFillDistanceVariable(DistanceVariable):
    
    def extendVariable(self, extensionVariable, order=2):
        if not hasattr(self, 'fineDistanceVariable'):
            self.fineDistanceVariable = DistanceVariable(mesh=self.mesh.fineMesh)
        if not hasattr(self, 'fineExtensionVariable'):
            self.fineExtensionVariable = CellVariable(mesh=self.mesh.fineMesh)
        self.fineDistanceVariable[:] = self(self.mesh.fineMesh.cellCenters)
        self.fineExtensionVariable[:] = extensionVariable(self.mesh.fineMesh.cellCenters)
        self.fineDistanceVariable.extendVariable(self.fineExtensionVariable, order=order)
        extensionVariable[:] = self.fineExtensionVariable(self.mesh.cellCenters)

    def calcDistanceFunction(self, order=2):
        if not hasattr(self, 'fineDistanceVariable'):
            self.fineDistanceVariable = DistanceVariable(mesh=self.mesh.fineMesh)
        self.fineDistanceVariable[:] = self(self.mesh.fineMesh.cellCenters)
        self.fineDistanceVariable.calcDistanceFunction(order=order)
        self[:] = self.fineDistanceVariable(self.mesh.cellCenters)

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test() 

