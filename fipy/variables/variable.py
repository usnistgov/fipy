"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "variable.py"
                                   created: 11/10/03 {3:15:38 PM} 
                               last update: 11/24/03 {10:29:55 PM} 
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
	self.name = name
	self.mesh = mesh
	self.array = Numeric.zeros([len(mesh.getCells())],'d')
        self.viewer = viewer
        if viewer != 'None':
            self.viewer.setVar(self)
	    
	self.setValue(value)

        if hasOld:
            self.old = Variable(name + "_old", mesh, value = self.getArray(), viewer = self.viewer, hasOld = 0)
        else:
            self.old = 'None'

    def plot(self):
        self.viewer.plot()
        
    def getMesh(self):
        return self.mesh

    def getArray(self):
        return self.array

    def getFaceArray(self):
        faces = self.mesh.getFaces()
        array = Numeric.zeros(len(faces),'d')
        for face in faces:
            array[face.getId()] = self.getFaceValue(face)
        return array

    def getGridArray(self):
        return self.mesh.makeGridData(self.array)
	
    def setValue(self,value,cells = ()):
	if cells == ():
            self.array[:] = value
        else:
            for cell in cells:
                self.array[cell.getId()] = value

    def getFaceValue(self,face):
        dAP = face.getCellDistance()
        dFP = face.getFaceToCellDistance()
        alpha = dFP / dAP
        id1 = face.getCellId(0)
        id2 = face.getCellId(1)
        return self.array[id1] * alpha + self.array[id2] * (1 - alpha)

    def getCellGradient(self,cell):
        grad = Numeric.zeros(self.mesh.getDim(),'d')
        for face in cell.getFaces():
            grad += self.getFaceValue(face) * face.getNormal(cell) * face.getArea()
        return grad/cell.getVolume()
        
    def getGradient(self):
        grad = Numeric.zeros((len(self.array),self.mesh.getDim()),'d')
        for cell in self.mesh.getCells():
            grad[cell.getId()] = self.getCellGradient(cell)
        return grad

    def getFaceGradient(self):
        faceGrad = Numeric.zeros((len(self.mesh.getFaces()),self.mesh.getDim()),'d')
        cellGrad = self.getGradient()
        for face in self.mesh.getFaces():
            dAP = face.getCellDistance()
            dFP = face.getFaceToCellDistance()
            alpha = dFP / dAP
            id1 = face.getCellId(0)
            id2 = face.getCellId(1)
            faceGrad[face.getId()] = cellGrad[id1] * alpha + cellGrad[id2] * (1 - alpha)
        return faceGrad

    def getCellGradientMagnitude(self,cell):
        grad  = self.getCellGradient(cell)
        mag = meshes.tools.sqrtDot(grad,grad)
        return mag
    
    def getGradientMagnitude(self):
        mag = Numeric.zeros(len(self.array),'d')
        for cell in self.mesh.getCells():
            mag[cell.getId()] = self.getCellGradientMagnitude(cell)
        return mag

    def getOld(self):
        if self.old == 'None':
            return self
        else:
            return self.old

    def updateOld(self):
        if self.old != 'None':
            self.old.setValue(self.array)

    
    
