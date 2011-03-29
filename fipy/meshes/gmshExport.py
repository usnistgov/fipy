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

class GmshExporter(object):

    def __init__(self, mesh, filename):
        self.mesh = mesh
        self.filename = filename
        self.outFile = open(self.filename, mode = 'w')

    def _getElementType(self, vertices, dimensions):
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

    def _orderVertices(self, vertexCoords, vertices):
        coordinates = numerix.take(vertexCoords, vertices, axis=1)
        centroid = numerix.add.reduce(coordinates, axis=1) / coordinates.shape[1]
        coordinates = coordinates - centroid[..., numerix.newaxis]

        # to prevent div by zero
        coordinates = numerix.where(coordinates == 0, 1.e-10, coordinates)

        # angles go from -pi / 2 to 3*pi / 2
        angles = numerix.arctan(coordinates[1] / coordinates[0]) \
                + numerix.where(coordinates[0] < 0, numerix.pi, 0) 
        sortorder = numerix.argsort(angles)
        return numerix.take(vertices, sortorder)
        

    def export(self):
        coords = self.mesh.vertexCoords
        dimensions, numNodes = coords.shape

        self._writeMeshFormat()
        self._writeNodes(coords, dimensions, numNodes)
        self._writeElements(coords, dimensions, numNodes)
        
        self.outFile.close()

    def _writeMeshFormat(self):
        versionNumber = 2.2
        sizeOfDouble = 8
        lines = ["$MeshFormat\n",
                 "%f 0 %d\n" % (versionNumber, sizeOfDouble),
                 "$EndMeshFormat\n"]
        self.outFile.writelines(lines)

    def _writeNodes(self, coords, dimensions, numNodes):
        self.outFile.write("$Nodes\n")
        self.outFile.write(str(numNodes) + '\n')

        for i in range(numNodes):
            self.outFile.write("%s %s %s " % (str(i + 1),
                                              str(coords[0, i]),
                                              str(coords[1, i])))
            if(dimensions == 2):
                self.outFile.write("0 \n")
            elif(dimensions == 3):
                self.outFile.write(str(coords[2, i]))
                self.outFile.write (" \n")
            else:
                raise MeshExportError, "Mesh has fewer than 2 or more than 3 dimensions" 

        self.outFile.write("$EndNodes\n")

    def _writeElements(self, coords, dimensions, numNodes):
        self.outFile.write("$Elements\n")

        faceVertexIDs = self.mesh.faceVertexIDs
        cellFaceIDs = self.mesh.cellFaceIDs
        numCells = cellFaceIDs.shape[1]
        self.outFile.write(str(numCells) + '\n')

        for i in range(numCells):
            ## build the vertex list
            vertexList = []
            for faceNum in cellFaceIDs[..., i]:
                """For more complicated meshes, some cells may have fewer
                faces than others. If this is the case, ignore the
                '--' entries."""
                if type(faceNum) == numerix.ma.core.MaskedConstant:
                    break
                for vertexNum in faceVertexIDs[..., faceNum]:
                    if vertexNum not in vertexList:
                        vertexList.append(vertexNum)

            if dimensions == 2:
                vertexList = self._orderVertices(coords, vertexList)

            numVertices = len(vertexList)
            elementType = self._getElementType(numVertices, dimensions)
            self.outFile.write("%s %s 0" % (str(i + 1), str(elementType)))

            for a in vertexList:
                self.outFile.write(" %s" % str(a + 1))
            self.outFile.write("\n")

        self.outFile.write("$EndElements\n") 

    def _test(self):
        """
        >>> import tempfile as tmp
        >>> from fipy.meshes.grid2D import Grid2D
        >>> g = Grid2D(nx = 10, ny = 10)
        >>> f = tmp.mktemp(".msh")
        >>> GmshExporter(g, f).export()

        >>> from fipy.meshes.uniformGrid2D import UniformGrid2D
        >>> ug = UniformGrid2D(nx = 10, ny = 10)
        >>> GmshExporter(ug, f).export()

        >>> from fipy.meshes import Tri2D
        >>> t = Tri2D(nx = 10, ny = 10)
        >>> GmshExporter(t, f).export()

        >>> concat = ug + (t + ([10], [0]))
        >>> GmshExporter(concat, f).export()

        >>> from fipy.meshes import Grid3D
        >>> g3d = Grid3D(nx=10, ny=10, nz=30)
        >>> GmshExporter(g3d, f).export()
        """

if __name__ == "__main__":
    from fipy.meshes import Grid2D
    from fipy.meshes import Tri2D
    from fipy.meshes import Grid3D
    from fipy.meshes import CylindricalGrid2D
    
    import tempfile as tmp
    import subprocess

    a = Grid2D(dx = 1.0, dy = 1.0, nx = 10, ny = 10)
    a2 = Grid2D(dx = 1.0, dy = 1.0, nx = 10, ny = 10) + ([10], [0])
    f1 = tmp.mktemp(".msh")
    GmshExporter(a, f1).export()

    b = Tri2D(dx = 1.0, dy = 1.0, nx = 10, ny = 10)
    f2 = tmp.mktemp(".msh")
    GmshExporter(b, f2).export()

    c = Grid3D(dx = 1.0, dy = 1.0, nx = 20, ny = 20, nz = 40)
    f3 = tmp.mktemp(".msh")
    GmshExporter(c, f3).export()
    subprocess.Popen(["gmsh", f3])
    raw_input("Grid3D... Press enter.")

    d = a + a2
    f4 = tmp.mktemp(".msh")
    GmshExporter(d, f4).export()
    subprocess.Popen(["gmsh", f4])
    raw_input("Concatenated grid... Press enter.")
 
    e = a + (b + ([0], [10]))
    f5 = tmp.mktemp(".msh")
    GmshExporter(e, f5).export()
    subprocess.Popen(["gmsh", f5])
    raw_input("Tri2D + Grid2D... Press enter.")

    cyl = CylindricalGrid2D(nx = 10, ny = 10)
    f6 = tmp.mktemp(".msh")
    GmshExporter(cyl, f6).export()
    subprocess.Popen(["gmsh", f6])
    raw_input("CylindricalGrid2D... Press enter.")

    
    
