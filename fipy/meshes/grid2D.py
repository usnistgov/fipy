#-----*-Pyth-*-
####################################################################
# PFM - Python-based phase field solver
#
# FILE: "grid2D.py"
#                                   created: 11/10/03 {3:30:42 PM} 
#                               last update: 11/13/03 {11:16:39 AM} 
# Author: Jonathan Guyer
# E-mail: guyer@nist.gov
#   mail: NIST
#    www: http://ctcms.nist.gov/
# 
#========================================================================
#This software was developed at the National Institute of Standards
#and Technology by employees of the Federal Government in the course
#of their official duties.  Pursuant to title 17 Section 105 of the
#United States Code this software is not subject to copyright
#protection and is in the public domain.  PFM is an experimental
#system.  NIST assumes no responsibility whatsoever for its use by
#other parties, and makes no guarantees, expressed or implied, about
#its quality, reliability, or any other characteristic.  We would
#appreciate acknowledgement if the software is used.
#
#This software can be redistributed and/or modified freely
#provided that any derivative works bear some notice that they are
#derived from it, and any modified versions bear some notice that
#they have been modified.
#========================================================================
# See the file "license.terms" for information on usage and  redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
# 
####################################################################
#----

import mesh
import vertex
import face
import cell

class Grid2D(mesh.Mesh):
    def __init__(self, dx, dy, nx, ny):
	self.vertices = self.createVertices(dx, dy, nx, ny)
	self.faceVertices,self.faceCells = self.createFaces(vertices, nx, ny)
	self.cells = self.createCells(faces, nx, ny)
	
	mesh.Mesh.__init__(self)
		
    def createVertices(self, dx, dy, nx, ny):
	vertices = Numeric.zeroes([(nx+1) * (ny+1),2],'d')
	for j in range(ny+1):
	    for	i in range(nx+1):
		vertices[i + (j * (nx+1)),0] = i * dx
		vertices[i + (j * (nx+1)),1] = j * dy
	return vertices
		    
    def createFaces(self, vertices, nx, ny):
	numFaces = (nx + 1) * ny + (ny + 1) * nx
	faceVertices = Numeric.zeroes([numFaces,2])
	faceCells = Numeric.zeroes([numFaces,2])

	for j in range(ny+1):
	    for i in range(nx):
		faceVertices[i +  (j * nx),0] = i + (j * nx)
		faceVertices[i +  (j * nx),1] = i + 1 + (j * nx)
		
	for j in range(ny):
	    for i in range(nx):
		faceCells[i +  (j * nx),0] = i + (j * nx)
		faceCells[i +  ((j+1) * nx),1] = i + (j * nx)
		
	for i in range(nx):
	    faceCells[i,0] =
	    faceCells[i + (ny+1) * nx,1] = 
		
	numHorzFaces = nx * (ny + 1)
	for j in range(ny):
	    for i in range(nx+1):
		faceVertices[numHorzFaces + (i * ny) + j,0] = (i * ny) + j
		faceVertices[numHorzFaces + (i * ny) + j,1] = (i * ny) + j + 1
		
	return (faceVertices,faceCells)
	
    def createCells(self, faces, nx, ny):
	cells = ()
	for j in range(ny):
	    for i in range(nx):
		cells += (cell.Cell((faces[i + j * nx],
			faces[i + (j+1) * nx],
			faces[nx * (ny + 1) + i + j * (nx + 1)],
			faces[nx * (ny + 1) + i + 1 + j * (nx + 1)])),)
	return cells
	
    def facesCellIDs(self):
	self.facesCellIDs = Numeric.zeroes([len(self.faces),2])
	for i in range(len(self.faces)):
	    self.facesCellIDs[i,0] = self.faces[i].cells()[0].id()
	    self.facesCellIDs[i,1] = self.faces[i].cells()[1].id()
	
    def cellsFaceIDs(self):
	pass
	
	
		
