#!/usr/bin/env python

##
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "gmshImport.py"
 #                                    created: 11/10/03 {2:44:42 PM}
 #                                last update: 10/22/04 {4:21:35 PM}
 #  Author: Alexander Mont <alexander.mont@nist.gov>
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2004-8-12 ADM 1.0 original
 # ###################################################################
 #


"""

This module takes a Gmsh output file (`.msh`) and converts it into a
FiPy mesh. This currently supports triangular and tetrahedral meshes
only.

Gmsh generates unstructured meshes, which may contain a significant
amount of non-orthogonality and it is very difficult to directly
control the amount of non-orthogonality simply by manipulating Gmsh
parameters. Therefore, it is necessary to take into account the
possibility of errors arising due to the non-orthogonality of the
mesh. To test the degree of error, an experiment was conducted using a
simple 1D diffusion problem with constant diffusion coefficients and
boundary conditions as follows: fixed value of 0 on the left side,
fixed value of 1 on the right side, and a fixed flux of 0 on the top
and bottom sides. The analytical solution is clearly a uniform
gradient going from left to right. this problem was implemented using
a Cartesian Grid2D mesh with each interior vertex displaced a short
distance in a random direction to create non-orthogonality. Then the
root-mean-square error was plotted against the root-mean-square
non-orthogonality. The error in each cell was calculated by simply
subtracting the analytical solution at each cell center from the
calculated value for that cell. The non-orthogonality in each cell is
the average, weighted by face area, of the sines of the angles between
the face normals and the line segments joining the cells. Thus, the
non-orthogonality of a cell can range from 0 (every face is orthogonal
to its corresponding cell-to-cell line segment) to 1 (only possible in
a degenerate case). This test was run using 500 separate 20x20 meshes
and 500 separate 10x10 meshes, each with the interior vertices moved
different amounts so as to created different levels of
non-orthogonality. The results are shown below.

Results for 20x20 mesh:

.. image:: images/orthoerrorgraph.pdf
   :height: 100
   :width: 200

Results for 10x10 mesh:

.. image:: images/orthoerrorcoarsegraph.pdf

It is clear from the graphs that finer meshes decrease the error due
to non-orthogonality, and that even with a reasonably coarse mesh the
error is quite low. However, note that this test is only for a simple
1D diffusion problem with a constant diffusion coefficient, and it is
unknown whether the results will be significantly different with more
complicated problems.

Test cases:

   >>> newmesh = GmshImporter3D('fipy/meshes/numMesh/testgmsh.msh')
   >>> print newmesh.getVertexCoords().tolist()
   [[0.0, 0.0, 0.0], [0.5, 0.5, 1.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0], [0.5, 0.5, 0.5]]

   >>> print newmesh.faceVertexIDs.tolist()
   [[2, 1, 0], [4, 1, 0], [4, 2, 0], [4, 2, 1], [3, 1, 0], [4, 3, 0], [4, 3, 1], [3, 2, 0], [4, 3, 2], [3, 2, 1]]

   >>> print newmesh.cellFaceIDs.tolist()
   [[0, 1, 2, 3], [4, 1, 5, 6], [7, 2, 5, 8], [9, 3, 6, 8]]

   >>> twomesh = GmshImporter2D('fipy/meshes/numMesh/GmshTest2D.msh')
   >>> print twomesh.getVertexCoords().tolist()
   [[0.0, 0.0], [1.0, 0.0], [0.5, 0.5], [0.0, 1.0], [1.0, 1.0], [0.5, 1.5], [0.0, 2.0], [1.0, 2.0]]
   
   >>> print twomesh.faceVertexIDs.tolist()
   [[2, 0], [0, 1], [1, 2], [0, 3], [3, 2], [1, 4], [4, 2], [4, 3], [3, 5], [5, 4], [3, 6], [6, 5], [5, 7], [7, 4], [7, 6]]
   
   >>> print twomesh.cellFaceIDs.tolist()
   [[0, 1, 2], [0, 3, 4], [2, 5, 6], [7, 4, 6], [7, 8, 9], [8, 10, 11], [12, 13, 9], [14, 11, 12]]
   
"""

__docformat__ = 'restructuredtext'

import Numeric
import MA
import mesh
import mesh2D
## from fipy.tools.profiler.profiler import Profiler
## from fipy.tools.profiler.profiler import calibrate_profiler

MeshImportError = "MeshImportError"

class DataGetter:

    def getData(self, fileName, dimensions):

        self.dimensions = dimensions
        self.fileName = fileName
        self.computeVertexCoords()
        self.computeCellVertexIDs()
        self.computeBaseFaceVertexIDs()
        self.computeFaceVertexIDs()
        self.computeCellFaceIDs()
        self.inFile.close()
        return (self.vertexCoords, self.faceVertexIDs, self.cellFaceIDs)

    def computeVertexCoords(self):

        dimensions = self.dimensions
        if (dimensions != 2 and dimensions != 3):
            raise MeshImportError, "Number of dimensions must be 2 or 3"
        
    ## initialize the file input stream
        inFile = open(self.fileName)
        a = inFile.readline() ## skip the $NOD

    ## get the vertex coordinates
        nodeToVertexIDdict = {}
        numVertices = int(inFile.readline())
        vertexCoords = Numeric.zeros((numVertices, dimensions))
        vertexCoords = vertexCoords.astype(Numeric.Float)
        for i in range(numVertices):
            currLine = inFile.readline()
            currLineArray = currLine.split()
            currLineArray[0] = int(currLineArray[0])
            for j in range(1, (dimensions + 1)):
                currLineArray[j] = float(currLineArray[j])
                nodeToVertexIDdict[currLineArray[0]] = i
            vertexCoords[i] = currLineArray[1:(dimensions + 1)]

        self.vertexCoords = vertexCoords
        self.inFile = inFile
        maxNode = max(nodeToVertexIDdict.keys())
        nodeToVertexIDs = Numeric.zeros((maxNode + 1,))
        for i in nodeToVertexIDdict.keys():
            nodeToVertexIDs[i] = nodeToVertexIDdict[i]
        self.nodeToVertexIDs = nodeToVertexIDs

    def computeCellVertexIDs(self):

        dimensions = self.dimensions
        inFile = self.inFile
        vertexCoords = self.vertexCoords
    
    ## get the elements
    ## note: all we care about are the three-dimensional elements (cells).
    ## note: so far this only supports tetrahedral and triangular meshes.
        a = inFile.readline() ## skip the $ENDNOD
        a = inFile.readline() ## skip the $ELM
        numElements = int(inFile.readline())
        numCells = 0
        maxLength = (6 + dimensions)
        elementArray = Numeric.zeros((numElements, maxLength))
        for i in range(numElements):
            currLineArrayInt = [int(x) for x in inFile.readline().split()]
            elementArray[i, :len(currLineArrayInt)] = currLineArrayInt
        validElementArray = Numeric.compress(elementArray[:, 1] == ((2 * dimensions) - 2), elementArray, 0)
        cellNodeIDs = validElementArray[:, 5:]
        cellVertexIDs = Numeric.take(self.nodeToVertexIDs, cellNodeIDs)        
        self.cellVertexIDs = cellVertexIDs
        self.numCells = len(cellVertexIDs)


    def computeBaseFaceVertexIDs(self):
        
        dimensions = self.dimensions
        cellVertexIDs = self.cellVertexIDs
    ## compute the face vertex IDs.
        cellFaceVertexIDs = Numeric.ones((self.numCells, dimensions + 1, dimensions))
        cellFaceVertexIDs = -1 * cellFaceVertexIDs

        if (dimensions == 3):
            cellFaceVertexIDs[:, 0, :] = cellVertexIDs[:, :3]
            cellFaceVertexIDs[:, 1, :] = Numeric.concatenate((cellVertexIDs[:, :2], cellVertexIDs[:, 3:]), axis = 1)
            cellFaceVertexIDs[:, 2, :] = Numeric.concatenate((cellVertexIDs[:, :1], cellVertexIDs[:, 2:]), axis = 1)
            cellFaceVertexIDs[:, 3, :] = cellVertexIDs[:, 1:]
        if (dimensions == 2):
            cellFaceVertexIDs[:, 0, :] = cellVertexIDs[:, :2]
##            cellFaceVertexIDs[:, 1, :] = Numeric.concatenate((cellVertexIDs[:, :1], cellVertexIDs[:, 2:]), axis = 1)
            cellFaceVertexIDs[:, 1, :] = Numeric.concatenate((cellVertexIDs[:, 2:], cellVertexIDs[:, :1]), axis = 1)
            cellFaceVertexIDs[:, 2, :] = cellVertexIDs[:, 1:]

        cellFaceVertexIDs = cellFaceVertexIDs[:, :, ::-1]
        self.unsortedBaseIDs = Numeric.reshape(cellFaceVertexIDs, (self.numCells * (dimensions + 1), dimensions))

        cellFaceVertexIDs = Numeric.sort(cellFaceVertexIDs, axis = 2)
        baseFaceVertexIDs = Numeric.reshape(cellFaceVertexIDs, (self.numCells * (dimensions + 1), dimensions))

        self.baseFaceVertexIDs = baseFaceVertexIDs       
        self.cellFaceVertexIDs = cellFaceVertexIDs

    def computeFaceVertexIDs(self):

        dimensions = self.dimensions
        faceStrToFaceIDs = {}
        faceStrToFaceIDsUnsorted = {}

        currIndex = 0

        for i in range(len(self.baseFaceVertexIDs)):
            listI = self.baseFaceVertexIDs[i]
            listJ = self.unsortedBaseIDs[i]
##            for i in self.baseFaceVertexIDs:
            if(not (faceStrToFaceIDs.has_key(listToString(listI)))):
                faceStrToFaceIDs[listToString(listI)] = currIndex
                faceStrToFaceIDsUnsorted[listToString(listJ)] = currIndex

                currIndex = currIndex + 1
        numFaces = currIndex
        faceVertexIDs = Numeric.zeros((numFaces, dimensions))
        for i in faceStrToFaceIDsUnsorted.keys():
            faceVertexIDs[faceStrToFaceIDsUnsorted[i], :] = stringToList(i)

        self.faceVertexIDs = faceVertexIDs
        self.faceStrToFaceIDs = faceStrToFaceIDs

    def computeCellFaceIDs(self):

        cellFaceIDs = Numeric.zeros(self.cellFaceVertexIDs.shape[:2])
        for i in range(len(self.cellFaceVertexIDs)):
            cell = self.cellFaceVertexIDs[i]
            for j in range(len(cell)):
                cellFaceIDs[i, j] = self.faceStrToFaceIDs[listToString(self.cellFaceVertexIDs[i, j])]
        self.cellFaceIDs = cellFaceIDs
            
class GmshImporter2D(mesh2D.Mesh2D):

    def __init__(self, filename):
        dg = DataGetter()
        a = dg.getData(filename, dimensions = 2)
        mesh2D.Mesh2D.__init__(self, a[0], a[1], a[2])

class GmshImporter3D(mesh.Mesh):

    def __init__(self, filename):
        dg = DataGetter()
        a = dg.getData(filename, dimensions = 3)
        mesh.Mesh.__init__(self, a[0], a[1], a[2])

def listToString(list):
    res = str(list[0])
    for i in list[1:]:
        res = res + ' ' + str(i)
    return res

def stringToList(string):
    return [int(x) for x in string.split(' ')]
    
def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    ##fudge = calibrate_profiler(10000)
    ##profile = Profiler('profile', fudge=fudge)
    ##newmesh = GmshImporter3D('untitled.msh')
    ##profile.stop()

    newmesh = GmshImporter3D('fipy/meshes/numMesh/testgmsh.msh')
    
    _test()
