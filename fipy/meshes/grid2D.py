#-----*-Pyth-*-
####################################################################
# PFM - Python-based phase field solver
#
# FILE: "grid2D.py"
#                                   created: 11/10/03 {3:30:42 PM} 
#                               last update: 11/10/03 {4:59:03 PM} 
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
	vertices = self.createVertices(dx, dy, nx, ny)
	faces = self.createFaces(vertices, nx, ny)
	cells = self.createCells(faces, nx, ny)
	
	mesh.Mesh.__init__(self, cells, faces, vertices)
		
    def createVertices(self, dx, dy, nx, ny):
	    vertices = ()
	    for j in range(ny+1):
		    for	i in range(nx+1):
			vertices += (vertex.Vertex((i * dx, j * dy)),)
	    return vertices
		
    def createFaces(self, vertices, nx, ny):
	faces = ()
	for j in range(ny+1):
	    for i in range(nx):
		faces += (face.Face((vertices[i + j * nx],vertices[i + 1 + j * nx])),)
	for j in range(ny):
	    for i in range(nx+1):
		faces += (face.Face((vertices[i * ny + j],vertices[i * ny + j + 1])),)
	return faces
	
    def createCells(self, faces, nx, ny):
	cells = ()
	for j in range(ny):
		for i in range(nx):
		    cells += (cell.Cell((faces[i + j * nx],
					    faces[i + (j+1) * nx],
					    faces[nx * (ny + 1) + i + j * (nx + 1)],
					    faces[nx * (ny + 1) + i + 1 + j * (nx + 1)])),)
	return cells
	
	
		
