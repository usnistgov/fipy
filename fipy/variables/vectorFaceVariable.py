#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "vectorFaceVariable.py"
 #                                    created: 12/9/03 {3:22:07 PM} 
 #                                last update: 4/4/05 {11:04:24 PM} 
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

import Numeric

from fipy.variables.variable import Variable

class VectorFaceVariable(Variable):
    def __init__(self,mesh,name = '',value=0., unit = None):

	array = Numeric.zeros([mesh._getNumberOfFaces(), mesh.getDim()],'d')
	
	Variable.__init__(self, mesh = mesh, name = name, value = value, unit = unit, array = array)

        self.indexAsFaceVar = {}

    def __call__(self, point = None, order = 0):
        if point != None:
            return self[self.getMesh()._getNearestCellID(point)]
        else:
            return Variable.__call__(self)
            
    def _getVariableClass(self):
	return VectorFaceVariable

    def getIndexAsFaceVariable(self, index):

        if not self.indexAsFaceVar.has_key(index):

            from faceVariable import FaceVariable
        
            class ItemAsVariable(FaceVariable):
                def __init__(self, var, index):
                    FaceVariable.__init__(self, mesh = var.getMesh())
                    self.var = self._requires(var)

                def _calcValue(self):
                    self.value = Numeric.array(self.var[:, index])

            return ItemAsVariable(self, index)

        else:
            
            return self.indexAsFaceVar[index]
            
    def getDivergence(self):
        if not hasattr(self, 'divergence'):
            from fipy.variables.addOverFacesVariable import AddOverFacesVariable
            self.divergence = AddOverFacesVariable(self.dot(self.getMesh()._getOrientedAreaProjections()))
            
        return self.divergence
        


	
