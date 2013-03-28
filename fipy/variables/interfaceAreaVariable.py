#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "interfaceAreaVariable.py"
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
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable
from fipy.tools.numerix import MA
from fipy.tools import numerix

class _InterfaceAreaVariable(CellVariable):
    def __init__(self, distanceVar):
        """
        Creates an `_InterfaceAreaVariable` object.

        :Parameters:
          - `distanceVar` : A `DistanceVariable` object.

        """
        CellVariable.__init__(self, distanceVar.mesh, hasOld=False)
        self.distanceVar = self._requires(distanceVar)

    def _calcValue(self):
        normals = numerix.array(MA.filled(self.distanceVar._cellInterfaceNormals, 0))
        areas = numerix.array(MA.filled(self.mesh._cellAreaProjections, 0))
        return numerix.sum(abs(numerix.dot(normals, areas)), axis=0)


    



            
            
        
                
