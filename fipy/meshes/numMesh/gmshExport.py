#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "gmshExport.py"
 #
 #  Author: Alexander Mont <alexander.mont@nist.gov>
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  gmshExport.py
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

"""
This module takes a FiPy mesh and creates a mesh file that can be opened in Gmsh.
"""

from fipy.tools import numerix
## from fipy.tools.profiler.profiler import Profiler
## from fipy.tools.profiler.profiler import calibrate_profiler

class MeshExportError(Exception):
    pass

def _getElementType(vertices, dimensions):
    if(vertices == 3 and dimensions == 2):
        return 2 ## triangle
    elif(vertices == 4 and dimensions == 2):
        return 3 ## quadrangle
    elif(vertices == 4 and dimensions == 3):
        return 4 ## tetrahedron
    elif(vertices == 8 and dimensions == 3):
        return 5 ## hexahedron
    elif(vertices == 6 and dimensions == 3):
        return 6 ## prism
    elif(vertices == 5 and dimensions == 3):
        return 7 ## pyramid
    else:
        raise MeshExportError, "Element type unsupported by Gmsh"

def _orderVertices(vertexCoords, vertices):
    coordinates = numerix.take(vertexCoords, vertices, axis=-1)
    centroid = numerix.add.reduce(coordinates) / coordinates.shape[0]
    coordinates = coordinates - centroid
    coordinates = numerix.where(coordinates == 0, 1.e-10, coordinates) ## to prevent division by zero
    angles = numerix.arctan(coordinates[:, 1] / coordinates[:, 0]) + numerix.where(coordinates[:, 0] < 0, numerix.pi, 0) ## angles go from -pi / 2 to 3*pi / 2
    sortorder = numerix.argsort(angles)
    return numerix.take(vertices, sortorder, axis=-1)
    

def exportAsMesh(mesh, filename):
    outFile = open(filename, mode = 'w')
    ## do the nodes
    outFile.write("$NOD\n")
    numNodes = mesh.getVertexCoords().shape[0]
    dimensions = mesh.getVertexCoords().shape[1]
    outFile.write(str(numNodes))
    outFile.write('\n')
    for i in range(numNodes):
        outFile.write(str(i + 1))
        outFile.write(' ')
        outFile.write(str(mesh.getVertexCoords()[i, 0]))
        outFile.write(' ') 
        outFile.write(str(mesh.getVertexCoords()[i, 1]))
        outFile.write(' ')
        if(dimensions == 2):
            outFile.write("0 \n")
        elif(dimensions == 3):
            outFile.write(str(mesh.getVertexCoords()[i, 2]))
            outFile.write (" \n")
        else:
            raise MeshExportError, "Mesh has fewer than 2 or more than 3 dimensions"
    outFile.write("$ENDNOD\n$ELM\n")
    ## do the elements
    faceVertexIDs = mesh.faceVertexIDs
    cellFaceIDs = mesh.cellFaceIDs
    numCells = cellFaceIDs.shape[0]
    outFile.write(str(numCells))
    outFile.write('\n')
    for i in range(numCells):
        ## build the vertex list
        vertexList = []
        for faceNum in cellFaceIDs[i]:
            for vertexNum in faceVertexIDs[faceNum]:
                if vertexNum not in vertexList:
                    vertexList = vertexList + [vertexNum]
        if(dimensions == 2):
            vertexList = _orderVertices(mesh.getVertexCoords(), vertexList)
        numVertices = len(vertexList)
        elementType = _getElementType(numVertices, dimensions)
        outFile.write(str(i + 1))
        outFile.write(' ')
        outFile.write(str(elementType))
        outFile.write(" 1 1 ")
        outFile.write(str(numVertices))
        for a in vertexList:
            outFile.write(' ')
            outFile.write(str(a + 1))
        outFile.write("\n")
    outFile.write("$ENDNOD\n")
    outFile.close()

if __name__ == "__main__":
    from fipy.meshes.grid2D import Grid2D
    from fipy.meshes.tri2D import Tri2D
    from fipy.meshes.grid3D import Grid3D
    import os
    ##a = Grid2D(dx = 1.0, dy = 1.0, nx = 10, ny = 10)
    ##exportAsMesh(a, "temp.msh")
    ##os.system("gmsh temp.msh &")
    ##b = Tri2D(dx = 1.0, dy = 1.0, nx = 10, ny = 10)
    ##exportAsMesh(b, "temp2.msh")
    ##os.system("gmsh temp2.msh &")
    fudge = calibrate_profiler(10000)
    profile = Profiler('profile', fudge=fudge)
    c = Grid3D(dx = 1.0, dy = 1.0, nx = 20, ny = 20, nz = 20)
    exportAsMesh(c, "temp3.msh")
    profile.stop()
    os.system("gmsh -v 0 temp3.msh &")
    
        
    
    
