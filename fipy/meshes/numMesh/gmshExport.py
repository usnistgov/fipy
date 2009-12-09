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
    coordinates = numerix.take(vertexCoords, vertices, axis=1)
    centroid = numerix.add.reduce(coordinates, axis=1) / coordinates.shape[1]
    coordinates = coordinates - centroid[..., numerix.newaxis]
    coordinates = numerix.where(coordinates == 0, 1.e-10, coordinates) ## to prevent division by zero
    angles = numerix.arctan(coordinates[1] / coordinates[0]) + numerix.where(coordinates[0] < 0, numerix.pi, 0) ## angles go from -pi / 2 to 3*pi / 2
    sortorder = numerix.argsort(angles)
    return numerix.take(vertices, sortorder)
    

def exportAsMesh(mesh, filename):
    outFile = open(filename, mode = 'w')
    ## do the nodes
    outFile.write("$NOD\n")
    coords = mesh.getVertexCoords()
    dimensions, numNodes = coords.shape
    outFile.write(str(numNodes))
    outFile.write('\n')
    for i in range(numNodes):
        outFile.write(str(i + 1))
        outFile.write(' ')
        outFile.write(str(coords[0, i]))
        outFile.write(' ') 
        outFile.write(str(coords[1, i]))
        outFile.write(' ')
        if(dimensions == 2):
            outFile.write("0 \n")
        elif(dimensions == 3):
            outFile.write(str(coords[2, i]))
            outFile.write (" \n")
        else:
            raise MeshExportError, "Mesh has fewer than 2 or more than 3 dimensions"
    outFile.write("$ENDNOD\n$ELM\n")
    ## do the elements
    faceVertexIDs = mesh._getFaceVertexIDs()
    cellFaceIDs = mesh._getCellFaceIDs()
    numCells = cellFaceIDs.shape[1]
    outFile.write(str(numCells))
    outFile.write('\n')
    for i in range(numCells):
        ## build the vertex list
        vertexList = []
        for faceNum in cellFaceIDs[..., i]:
            for vertexNum in faceVertexIDs[..., faceNum]:
                if vertexNum not in vertexList:
                    vertexList = vertexList + [vertexNum]
        if(dimensions == 2):
            vertexList = _orderVertices(coords, vertexList)
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
    
    import tempfile
    import subprocess

    a = Grid2D(dx = 1.0, dy = 1.0, nx = 10, ny = 10)
    f = tempfile.NamedTemporaryFile(suffix=".msh")
    f.close()
    exportAsMesh(a, f.name)
    subprocess.Popen(["gmsh", f.name])
    b = Tri2D(dx = 1.0, dy = 1.0, nx = 10, ny = 10)
    f = tempfile.NamedTemporaryFile(suffix=".msh")
    f.close()
    exportAsMesh(b, f.name)
    subprocess.Popen(["gmsh", f.name])
#     fudge = calibrate_profiler(10000)
#     profile = Profiler('profile', fudge=fudge)
    c = Grid3D(dx = 1.0, dy = 1.0, nx = 20, ny = 20, nz = 20)
    f = tempfile.NamedTemporaryFile(suffix=".msh")
    f.close()
    exportAsMesh(c, f.name)
#     profile.stop()
    subprocess.Popen(["gmsh", "-v", "0", f.name])
    
        
    
    
