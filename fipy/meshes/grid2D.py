"""
#-----*-Pyth-*-
####################################################################
# PFM - Python-based phase field solver
#
# FILE: "grid2D.py"
#                                   created: 11/10/03 {3:30:42 PM} 
#                               last update: 11/13/03 {11:22:42 AM} 
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
"""
import mesh
import vertex
import face
import cell

class Grid2D(mesh.Mesh):
    """
    NUMBERING SYSTEM
    ----------------

    CELLS
    -----
    nx=5
    ny=3
    
    *************************************
    *      *      *       *      *      *
    * 10   * 11   * 12    * 13   * 14   *
    *************************************
    *      *      *       *      *      *
    * 5    * 6    * 7     * 8    * 9    *
    *************************************
    *      *      *       *      *      *
    * 0    * 1    * 2     * 3    * 4    *
    *************************************

    FACES
    -----
    ***15******16*****17******18****19***
    *      *      *       *      *      *
    32    33      34      35     36     37
    ***10******11*****12******13*****14**
    *      *      *       *      *      *
    26     27     28      29     30     31
    ***5*******6******7*******8******9***
    *      *      *       *      *      *
    20     21     22      23     24     25
    ***0*******1******2*******3******4***

    VERTICES
    --------

    18*****19*****20******21*****22****23
    *      *      *       *      *      *
    *      *      *       *      *      *
    12*****13*****14******15*****16****17
    *      *      *       *      *      *
    *      *      *       *      *      *
    6******7******8*******9******10****11
    *      *      *       *      *      *
    *      *      *       *      *      *
    0******1******2*******3******4******5
    """
    
    def __init__(self, dx, dy, nx, ny):
        self.nx=nx
        self.ny=ny
        self.dx=dx
        self.dy=dy
	vertices = self.createVertices()
	faces = self.createFaces(vertices)
	cells = self.createCells(faces)
        interiorFaces = self.createInteriorFaces(faces)
	
	mesh.Mesh.__init__(self, cells, faces, interiorFaces, vertices)
		
    def createVertices(self):
	vertices = ()
        ny=self.ny
        nx=self.nx
        dx=self.dx
        dy=self.dy
	for j in range(ny+1):
	    for	i in range(nx+1):
		vertices += (vertex.Vertex((i * dx, j * dy)),)
        return vertices	
		    
    def createFaces(self, vertices):
        nx=self.nx
        ny=self.ny
	faces = ()
	for j in range(ny+1):
	    for i in range(nx):
		faces += (face.Face((vertices[i + j * nx],vertices[i + 1 + j * nx])),)
	for j in range(ny):
	    for i in range(nx+1):
		faces += (face.Face((vertices[i * ny + j],vertices[i * ny + j + 1])),)
	return faces
	
    def createCells(self, faces):
        nx=self.nx
        ny=self.ny
	cells = ()
	for j in range(ny):
	    for i in range(nx):
                id = j * nx + i
		cells += (
                    cell.Cell(
                    (faces[i + j * nx],
                    faces[i + (j+1) * nx],
                    faces[nx * (ny + 1) + i + j * (nx + 1)],
                    faces[nx * (ny + 1) + i + 1 + j * (nx + 1)]),id
                    ),
                    )                
	return cells

    def createInteriorFaces(self,faces):
        interiorFaces = ()
        for face in faces:
            if len(face.getCells()) == 1:
                interiorFaces += (face,)
        return interiorFaces

    def getFacesLeft(self):
        nx=self.nx
        ny=self.ny
        facesLeft = ()
        for i in range(ny):
            facesLeft += (self.faces[nx * (ny + 1) + (nx + 1) * i],)
        return facesLeft
        
    def getFacesRight(self):
        nx=self.nx
        ny=self.ny
        facesRight = ()
        for i in range( ny):
            facesRight += (self.faces[nx * (ny + 1) + (nx + 1) * (i + 1) - 1],)
        return facesRight

    def getFacesTop(self):
        nx=self.nx
        ny=self.ny
        return self.faces[nx*ny:nx*ny+nx]

    def getFacesBottom(self):
        nx=self.nx
        return self.faces[0:nx]
                       
        
        
