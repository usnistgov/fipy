"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "variable.py"
                                   created: 11/10/03 {3:15:38 PM} 
                               last update: 11/30/03 {12:56:28 AM} 
 Author: Jonathan Guyer
 E-mail: guyer@nist.gov
 Author: Daniel Wheeler
 E-mail: daniel.wheeler@nist.gov
   mail: NIST
    www: http://ctcms.nist.gov
 
========================================================================
This software was developed at the National Institute of Standards
and Technology by employees of the Federal Government in the course
of their official duties.  Pursuant to title 17 Section 105 of the
United States Code this software is not subject to copyright
protection and is in the public domain.  PFM is an experimental
system.  NIST assumes no responsibility whatsoever for its use by
other parties, and makes no guarantees, expressed or implied, about
its quality, reliability, or any other characteristic.  We would
appreciate acknowledgement if the software is used.

This software can be redistributed and/or modified freely
provided that any derivative works bear some notice that they are
derived from it, and any modified versions bear some notice that
they have been modified.
========================================================================
 
 Description: 

 History

 modified   by  rev reason
 ---------- --- --- -----------
 2003-11-10 JEG 1.0 original
###################################################################
"""

import meshes.tools
import Numeric

class Variable:
    
    def __init__(self, name, mesh, value=0.,viewer = 'None', hasOld = 1):
	self.fields = {}
	
	self.name = name
	self.mesh = mesh
	self.array = Numeric.zeros([len(mesh.getCells())],'d')
        if viewer != 'None':
	    self.viewer = viewer(var = self)
	else:
	    self.viewer = 'None'
	    
	self.setValue(value)

        if hasOld:
            self.old = Variable(name + "_old", mesh, value = self.getArray(), viewer = viewer, hasOld = 0)
        else:
            self.old = 'None'

    def plot(self):
	if self.viewer != 'None':
	    self.viewer.plot()
        
    def getMesh(self):
        return self.mesh

    def getArray(self):
        return self.array

    def getGridArray(self):
        return self.mesh.makeGridData(self.array)
	
    def setValue(self,value,cells = ()):
	if cells == ():
            self.array[:] = value
        else:
            for cell in cells:
                self.array[cell.getId()] = value
	self.refresh()

    def refresh(self):
	if self.fields.has_key('faceValues'):
	    self.calcFaceValues()
	if self.fields.has_key('gradient'):
	    self.calcGradient()
	
    def getFaceValue(self,face):
	return self.getFaceValues()[face.getId()]
	    
    def getFaceValues(self):
	if not self.fields.has_key('faceValues'):
	    self.calcFaceValues()
	return self.fields['faceValues']

    def calcFaceValues(self):
	dAP = self.mesh.getCellDistances()
	dFP = self.mesh.getFaceToCellDistances()
	alpha = dFP / dAP
	id1, id2 = self.mesh.getAdjacentCellIDs()
	self.fields['faceValues'] = Numeric.take(self.array, id1) * alpha + Numeric.take(self.array, id2) * (1 - alpha)
	    
    def calcFaceGradientContributions(self):
	areas = self.mesh.getAreaProjections()
	faceValues = Numeric.reshape(self.getFaceValues(), (len(areas),1)) 
	
	self.faceGradientContributions = areas * faceValues
    
    def getCellGradient(self,cell):
	contributions = Numeric.take(self.faceGradientContributions, cell.getFaceIDs())
	
	return Numeric.sum(cell.getFaceOrientations() * contributions) / cell.getVolume()
	
    def getGradient(self):
	if not self.fields.has_key('gradient'):
	    self.calcGradient()
	return self.fields['gradient']
	    
    def calcGradient(self):
	self.calcFaceGradientContributions()
        grad = Numeric.zeros((len(self.array),self.mesh.getDim()),'d')
        for cell in self.mesh.getCells():
            grad[cell.getId()] = self.getCellGradient(cell)
	self.fields['gradient'] = grad

    def getFaceGradient(self):
	cellGrad = self.getGradient()
	alpha = Numeric.zeros((len(self.mesh.getFaces()),1),'d')
	dAP = self.mesh.getCellDistances()
	dFP = self.mesh.getFaceToCellDistances()
	alpha[:,0] = dFP / dAP
	id1, id2 = self.mesh.getAdjacentCellIDs()
	return Numeric.take(cellGrad, id1) * alpha + Numeric.take(cellGrad, id2) * (1 - alpha)

    def getGradientMagnitude(self):
	grad = self.getGradient()
	return Numeric.sqrt(Numeric.sum(grad*grad,1))
	
    def getOld(self):
        if self.old == 'None':
            return self
        else:
            return self.old

    def updateOld(self):
        if self.old != 'None':
            self.old.setValue(self.array)

    
    
