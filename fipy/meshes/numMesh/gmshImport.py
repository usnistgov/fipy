#!/usr/bin/env python

##
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "gmshImport.py"
 #
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
 # ###################################################################
 #

r"""

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

.. image:: fipy/meshes/numMesh/orthoerrorgraph.*
   :height: 100px
   :width: 200px

Results for 10x10 mesh:

.. image:: fipy/meshes/numMesh/orthoerrorcoarsegraph.*
    :height: 100px
    :width: 200px

It is clear from the graphs that finer meshes decrease the error due
to non-orthogonality, and that even with a reasonably coarse mesh the
error is quite low. However, note that this test is only for a simple
1D diffusion problem with a constant diffusion coefficient, and it is
unknown whether the results will be significantly different with more
complicated problems.

Test cases:

   >>> newmesh = GmshImporter3D('fipy/meshes/numMesh/testgmsh.msh')
   >>> print newmesh.getVertexCoords()
   [[ 0.   0.5  1.   0.5  0.5]
    [ 0.   0.5  0.   1.   0.5]
    [ 0.   1.   0.   0.   0.5]]

   >>> print newmesh._getFaceVertexIDs()
   [[2 4 4 4 3 4 4 3 4 3]
    [1 1 2 2 1 3 3 2 3 2]
    [0 0 0 1 0 0 1 0 2 1]]

   >>> print newmesh._getCellFaceIDs()
   [[0 4 7 9]
    [1 1 2 3]
    [2 5 5 6]
    [3 6 8 8]]

   >>> mesh = GmshImporter2DIn3DSpace('fipy/meshes/numMesh/GmshTest2D.msh')
   >>> print mesh.getVertexCoords()
   [[ 0.   1.   0.5  0.   1.   0.5  0.   1. ]
    [ 0.   0.   0.5  1.   1.   1.5  2.   2. ]
    [ 0.   0.   0.   0.   0.   0.   0.   0. ]]

   >>> mesh = GmshImporter2D('fipy/meshes/numMesh/GmshTest2D.msh')
   >>> print mesh.getVertexCoords()
   [[ 0.   1.   0.5  0.   1.   0.5  0.   1. ]
    [ 0.   0.   0.5  1.   1.   1.5  2.   2. ]]

   >>> print mesh._getFaceVertexIDs()
   [[2 0 1 0 3 1 4 4 3 5 3 6 5 7 7]
    [0 1 2 3 2 4 2 3 5 4 6 5 7 4 6]]
   
   >>> print (mesh._getCellFaceIDs() == [[0, 0, 2, 7, 7, 8, 12, 14],
   ...                                   [1, 3, 5, 4, 8, 10, 13, 11],
   ...                                   [2, 4, 6, 6, 9, 11, 9, 12]]).flatten().all()
   True
   
The following test case is to test the handedness of the mesh to check
it does not return negative volumes. Firstly we set up a list with
tuples of strings to be read by gmsh. The list provide instuctions to
gmsh to form a circular mesh.

   >>> cellSize = 0.7
   >>> radius = 1.
   >>> lines = ['cellSize = ' + str(cellSize) + ';\n',
   ...           'radius = ' + str(radius) + ';\n',
   ...           'Point(1) = {0, 0, 0, cellSize};\n',
   ...           'Point(2) = {-radius, 0, 0, cellSize};\n',
   ...           'Point(3) = {0, radius, 0, cellSize};\n',
   ...           'Point(4) = {radius, 0, 0, cellSize};\n',
   ...           'Point(5) = {0, -radius, 0, cellSize};\n',
   ...           'Circle(6) = {2, 1, 3};\n',
   ...           'Circle(7) = {3, 1, 4};\n',
   ...           'Circle(8) = {4, 1, 5};\n',
   ...           'Circle(9) = {5, 1, 2};\n',
   ...           'Line Loop(10) = {6, 7, 8, 9} ;\n',
   ...           'Plane Surface(11) = {10};\n']

Check that the sign of the mesh volumes is correct

   >>> mesh = GmshImporter2D(lines)
   >>> print mesh.getCellVolumes()[0] > 0
   1
      
Reverse the handedness of the mesh and check the sign

   >>> lines[7:12] = ['Circle(6) = {3, 1, 2};\n',
   ...                'Circle(7) = {4, 1, 3};\n',
   ...                'Circle(8) = {5, 1, 4};\n',
   ...                'Circle(9) = {2, 1, 5};\n',
   ...                'Line Loop(10) = {9, 8, 7, 6};\n',]

   >>> mesh = GmshImporter2D(lines)
   >>> print mesh.getCellVolumes()[0] > 0
   1
   
"""

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools.numerix import MA
import mesh
import mesh2D
## from fipy.tools.profiler.profiler import Profiler
## from fipy.tools.profiler.profiler import calibrate_profiler
import os

class MeshImportError(Exception):
    pass

class _DataGetter:
    def __init__(self, filename, dimensions, coordDimensions = None):
        self.coordDimensions = coordDimensions or dimensions

        if (dimensions != 2 and dimensions != 3):
            raise MeshImportError, "Number of dimensions must be 2 or 3"
        self.dimensions = dimensions
        self.filename = filename
        
    def getData(self):
        self.inFile = open(self.filename)
        self.fileType = self.getFileType() #gets version of gmsh, I think
        #vertexCoords are x,y coords of nodes/vertices from gmsh file 
        vertexCoords = self._calcVertexCoords(self.coordDimensions)

        self._calcCellVertexIDs()
        self._calcBaseFaceVertexIDs()
        faceVertexIDs = self._calcFaceVertexIDs()#reads nodes/vertices from gmsh file 
        cellFaceIDs = self._calcCellFaceIDs()

        self.inFile.close()

        return {
            'vertexCoords': vertexCoords,
            'faceVertexIDs': faceVertexIDs,
            'cellFaceIDs': cellFaceIDs
            }
            
    def getFileType(self):
        data = self.getTagData("$MeshFormat", "$EndMeshFormat")
        if data is None:
            return 1.0
        else:
            return float(data[0].split()[0])
        
    def getTagData(self, begin, end):
        self.inFile.seek(0)
        
        for line in self.inFile:
            if begin in line:
                data = []
                for subline in self.inFile:
                    if end in subline:
                        return data
                    data.append(subline)
                raise EOFError, "No matching '%s' for '%s'" % (end, begin)
                
        return None

    def _calcVertexCoords(self, coordDimensions):
        if self.fileType == 1:
            nodeLines = self.getTagData("$NOD", "$ENDNOD")
        else:
            nodeLines = self.getTagData("$Nodes", "$EndNodes")

    ## get the vertex coordinates
        nodeToVertexIDdict = {}
        
        numVertices = int(nodeLines[0])
        if numVertices != len(nodeLines[1:]):
            raise IndexError, "Number of nodes (%d) does not match number promised (%d)" % (numVertices, len(nodeLines[1:]))

        vertexCoords = []
        for node, i in zip(nodeLines[1:], range(len(nodeLines[1:]))):
            nodeInfo = node.split()
            nodeToVertexIDdict[int(nodeInfo[0])] = i
            vertexCoords.append([float(n) for n in nodeInfo[1:]])

        vertexCoords = numerix.array(vertexCoords, 'd')
        
        maxNode = max(nodeToVertexIDdict.keys())
        nodeToVertexIDs = numerix.zeros((maxNode + 1,))
        for i in nodeToVertexIDdict.keys():
            nodeToVertexIDs[i] = nodeToVertexIDdict[i]
        self.nodeToVertexIDs = nodeToVertexIDs
        return vertexCoords[:,:coordDimensions].swapaxes(0,1)
        
    def _calcCellVertexIDs(self):
        """
        Get the elements.
        
        .. note:: all we care about are the three-dimensional elements (cells).
        
        .. note:: so far this only supports tetrahedral and triangular meshes.
        """
        if self.fileType == 1:
            elementLines = self.getTagData("$ELM", "$ENDELM")
        else:
            elementLines = self.getTagData("$Elements", "$EndElements")

        numElements = int(elementLines[0])
        if numElements != len(elementLines[1:]):
            raise IndexError, "Number of elements (%d) does not match number promised (%d)" % (numElements, len(elementLines[1:]))
                            
        cellNodeIDs = []
        for element in elementLines[1:]:
            elementInfo = [int(x) for x in element.split()]
            elementType = elementInfo[1]
            if elementType in (1, 15):
                continue
            elif elementType in (2, 4):
                if ((self.dimensions == 2 and elementType == 4) 
                    or (self.dimensions == 3 and elementType == 2)):
                    continue
                                  
                if self.fileType == 1:
                    numNodes = elementInfo[4]
                    skip = 5
                else:
                    if elementType == 2:
                        numNodes = 3
                    else:
                        numNodes = 4
                    skip = 3 + elementInfo[2] 
                    
                if len(elementInfo) != skip + numNodes:
                    raise IndexError, "Number of nodes (%d) not as expected (%d) for element type %d" % (len(elementInfo) - skip, numNodes, elementType)

                cellNodeIDs.append(elementInfo[skip:])
            elif elementType in (3,):
                if self.fileType == 2:
                    numNodes = 4
                    #skip is the number of columns to pass over to get to vertex info.
                    skip = 3 + elementInfo[2]
                else:
                    raise TypeError, "don't know how to handle quadralaterals in version 1. files"
                cellNodeIDs.append(elementInfo[skip:])
            else:
                raise TypeError, "Can't understand element type %d. Only triangle (2) or tetrahedron (4) are allowed" % elementInfo[1]
                
        self.cellVertexIDs = numerix.take(self.nodeToVertexIDs, 
                                          numerix.array(cellNodeIDs)).swapaxes(0,1)       
        self.numCells = self.cellVertexIDs.shape[-1]

    def _calcBaseFaceVertexIDs(self):
        
        cellVertexIDs = self.cellVertexIDs
    ## compute the face vertex IDs.
        ### this assumes triangular grid
        #cellFaceVertexIDs = numerix.ones((self.dimensions, self.dimensions + 1, self.numCells))
        cellFaceVertexIDs = numerix.ones((self.dimensions,len(cellVertexIDs), self.numCells))
        cellFaceVertexIDs = -1 * cellFaceVertexIDs

        if (self.dimensions == 3):
            cellFaceVertexIDs[:, 0, :] = cellVertexIDs[:3]
            cellFaceVertexIDs[:, 1, :] = numerix.concatenate((cellVertexIDs[:2], cellVertexIDs[3:]), axis = 0)
            cellFaceVertexIDs[:, 2, :] = numerix.concatenate((cellVertexIDs[:1], cellVertexIDs[2:]), axis = 0)
            cellFaceVertexIDs[:, 3, :] = cellVertexIDs[1:]
        elif (self.dimensions == 2):#define face with vertex pairs
            ###This isn't very general.
            ###Would be nice to allow cells with different number of faces. 
            if len(cellVertexIDs)==3:
                cellFaceVertexIDs[:, 0, :] = cellVertexIDs[:2]
                cellFaceVertexIDs[:, 1, :] = numerix.concatenate((cellVertexIDs[2:], cellVertexIDs[:1]), axis = 0)
                cellFaceVertexIDs[:, 2, :] = cellVertexIDs[1:]
            elif len(cellVertexIDs)==4:
                cellFaceVertexIDs[:, 0, :] = cellVertexIDs[0:2]
                cellFaceVertexIDs[:, 1, :] = cellVertexIDs[1:3]
                cellFaceVertexIDs[:, 2, :] = cellVertexIDs[2:4]
                cellFaceVertexIDs[:, 3, :] = numerix.concatenate((cellVertexIDs[3:], cellVertexIDs[:1]), axis = 0)

        cellFaceVertexIDs = cellFaceVertexIDs[::-1]#reverses order of vertex pair

        #self.unsortedBaseIDs = numerix.reshape(cellFaceVertexIDs.swapaxes(1,2), 
        #                                       (self.dimensions, 
        #                                        self.numCells * (self.dimensions + 1)))
        
        self.unsortedBaseIDs = numerix.reshape(cellFaceVertexIDs.swapaxes(1,2), 
                                               (self.dimensions, 
                                                self.numCells * (len(cellVertexIDs))))
        cellFaceVertexIDs = numerix.sort(cellFaceVertexIDs, axis=0)
        baseFaceVertexIDs = numerix.reshape(cellFaceVertexIDs.swapaxes(1,2), 
                                            (self.dimensions, 
                                             self.numCells * (len(cellVertexIDs))))

        self.baseFaceVertexIDs = baseFaceVertexIDs       
        self.cellFaceVertexIDs = cellFaceVertexIDs

    def _calcFaceVertexIDs(self):

        self.faceStrToFaceIDs = {}
        faceStrToFaceIDsUnsorted = {}

        currIndex = 0

        for i in range(self.baseFaceVertexIDs.shape[-1]):
            listI = self.baseFaceVertexIDs[:,i]
            listJ = self.unsortedBaseIDs[:,i]

            key = ' '.join([str(i) for i in listI])
            if(not (self.faceStrToFaceIDs.has_key(key))):
                self.faceStrToFaceIDs[key] = currIndex
                faceStrToFaceIDsUnsorted[' '.join([str(j) for j in listJ])] = currIndex

                currIndex = currIndex + 1
        numFaces = currIndex
        faceVertexIDs = numerix.zeros((self.dimensions, numFaces))
        for i in faceStrToFaceIDsUnsorted.keys():
            faceVertexIDs[:, faceStrToFaceIDsUnsorted[i]] = [int(x) for x in i.split(' ')]

        return faceVertexIDs

    def _calcCellFaceIDs(self):

        cellFaceIDs = numerix.zeros(self.cellFaceVertexIDs.shape[1:])
        for j in range(self.cellFaceVertexIDs.shape[-1]):
            cell = self.cellFaceVertexIDs[...,j]
            for i in range(cell.shape[-1]):
                cellFaceIDs[i, j] = self.faceStrToFaceIDs[' '.join([str(k) for k in self.cellFaceVertexIDs[:,i, j]])]
        return cellFaceIDs

class MshFile:
    def __init__(self, arg):

        if '.msh' in arg:
            self.mshfile = arg
            self.deletemshfile = False
        else:

            import tempfile
                   
            if not ('.geo' in arg or '.gmsh' in arg):
                (f, geofile) = tempfile.mkstemp('.geo')
                file = open(geofile, 'w')
                file.writelines(arg)
                file.close()
                os.close(f)
                self.deletegeofile = True
            else:
                self.deletegeofile = False
                geofile = arg
                
            (f, self.mshfile) = tempfile.mkstemp('.msh')

            os.system('gmsh ' + geofile + ' -2 -v 0 -format msh -o ' + self.mshfile)

            os.close(f)

            if self.deletegeofile:
                os.remove(geofile)

            self.deletemshfile = True

    def getFilename(self):
        return self.mshfile

    def remove(self):
        if self.deletemshfile:
            os.remove(self.mshfile)
                 
class GmshImporter2D(mesh2D.Mesh2D):

    def __init__(self, arg, coordDimensions=2):
        mshfile = MshFile(arg)
        mesh2D.Mesh2D.__init__(self, **_DataGetter(mshfile.getFilename(), dimensions=2, coordDimensions=coordDimensions).getData())
        mshfile.remove()
        
    def getCellVolumes(self):
        return abs(mesh2D.Mesh2D.getCellVolumes(self))

class GmshImporter2DIn3DSpace(GmshImporter2D):
    def __init__(self, arg):
        GmshImporter2D.__init__(self, arg, coordDimensions=3)

class GmshImporter3D(mesh.Mesh):
    """
        >>> mesh = GmshImporter3D('fipy/meshes/numMesh/testgmsh.msh')
    """

    def __init__(self, arg):
        mshfile = MshFile(arg)
        mesh.Mesh.__init__(self, **_DataGetter(mshfile.getFilename(), dimensions=3).getData())
        mshfile.remove()

    def getCellVolumes(self):
        return abs(mesh.Mesh.getCellVolumes(self))
    
def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    ##fudge = calibrate_profiler(10000)
    ##profile = Profiler('profile', fudge=fudge)
    ##newmesh = GmshImporter3D('untitled.msh')
    ##profile.stop()

    _test()
