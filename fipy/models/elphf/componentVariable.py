#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "componentVariable.py"
 #                                    created: 12/18/03 {12:18:05 AM} 
 #                                last update: 1/16/04 {11:00:24 AM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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

from variables.cellVariable import CellVariable
from tools.dimensions import physicalField

class ComponentVariable(CellVariable):
    def __init__(self, mesh, parameters, value=0., hasOld = 1):
	if parameters.has_key('name'):
	    name = parameters['name']
	else:
	    name = ''
# 	value = physicalField.PhysicalField(value)
	value = physicalField.Scale(value, "MOLARVOLUME**-1")
	CellVariable.__init__(self, mesh = mesh, name = name, value = value, hasOld = hasOld)
	self.parameters = parameters
## 	self.standardPotential = physicalField.PhysicalField(parameters['standard potential'])
## 	self.barrierHeight = physicalField.PhysicalField(parameters['barrier height'])
 	self.standardPotential = physicalField.Scale(parameters['standard potential'], "ENERGY")
 	self.barrierHeight = physicalField.Scale(parameters['barrier height'], "ENERGY")
	if self.parameters.has_key('valence'):
	    self.valence = self.parameters['valence']
	else:
	    self.valence = 0
	if self.parameters.has_key('diffusivity'):
 	    self.diffusivity = physicalField.Scale(parameters['diffusivity'], "LENGTH**2/TIME")
## 	    self.diffusivity = physicalField.PhysicalField(parameters['diffusivity'])
	else:
	    self.diffusivity = 0
	
    def getStandardPotential(self):
	return self.standardPotential
	
    def getBarrierHeight(self):
	return self.barrierHeight
	
    def getValence(self):
	return self.valence
	
    def getDiffusivity(self):
	return self.diffusivity