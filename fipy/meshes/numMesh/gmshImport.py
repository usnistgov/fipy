#!/usr/bin/env python
 
## -*-Pyth-*-
# ###################################################################
#  FiPy - a finite volume PDE solver in Python
#
#  FILE: "gmshImport.py"
#
#  Author: James O'Beirne <james.obeirne@nist.gov>
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
#  See the file "license.terms" for information on usage and
#  redistribution
#  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
# 
# ###################################################################
##

from fipy.tools import numerix as nx
import mesh
import mesh2D
import os
import tempfile

class MshFile:
    """
    Class responsible for parsing a Gmsh file and then readying
    its contents for use by a `Mesh` constructor.

    Does not support gmsh versions < 2.
    """
    def __init__(self, filename, dimensions, coordDimensions=None):
        """
        Isolates relevant data into two files, stores in 
        `self.nodesFile` for $Nodes,
        `self.elemsFile` for $Elements.

        :Parameters:
          - `filename`: a string indicating gmsh output file
          - `dimensions`: an integer indicating dimension of mesh
          - `coordDimension`: an integer indicating dimension of shapes
        """
        
        self.coordDimensions  = coordDimensions or dimensions
        self.dimensions       = dimensions
        self.filename         = self._parseFilename(filename)
        self.numFacesForShape = {2: 3, # triangle:   3 sides
                                 3: 4, # quadrangle: 4 sides
                                 4: 4} # tet:        4 sides

        f = open(self.filename, "r") # open the msh file

        # five lines of muck here allow most of the other methods to be
        # free of side-effects.
        self.version, self.fileType, self.dataSize = self._getMetaData(f)
        self.nodesFile = self._isolateData("Nodes", f)
        self.elemsFile = self._isolateData("Elements", f)

        self.vertexCoords, self.vertexMap = self._vertexCoordsAndMap()
        self.facesToV, self.cellsToF  = self._parseElements(self.vertexMap)

    def _parseFilename(self, fname, gmshFlags="-2 -v 0 -format msh"):
        """
        If we're being passed a .msh file, leave it be. Otherwise,
        we've gotta compile a .msh file from either (i) a .geo file, 
        or (ii) a gmsh script passed as a string.
        """
        lowerFname = fname.lower()
        if '.msh' in lowerFname:
            return fname
        else:
            if '.geo' in lowerFname or '.gmsh' in lowerFname:
                geoFile = fname
            else: # fname must be a full script, not a file
                (f, geoFile) = tempfile.mkstemp('.geo')
                file = open(geoFile, 'w')
                file.writelines(fname)
                file.close(); os.close(f)

            (f, mshFile) = tempfile.mkstemp('.msh')
            os.system('gmsh %s %s -o %s' \
                      % (geoFile, gmshFlags, mshFile))
            os.close(f)

            return mshFile
         
    def _getMetaData(self, f):
        """
        Extracts gmshVersion, file-type, and data-size in that
        order.
        """
        self._seekForHeader("MeshFormat", f)
        metaData = f.readline().split()
        f.seek(0)
        return [float(x) for x in metaData]

    def _isolateData(self, title, f):
        """
        Gets all data between $[title] and $End[title], writes
        it out to its own file.
        """
        newF = tempfile.TemporaryFile()

        # seek to section header
        self._seekForHeader(title, f)
        
        # extract the actual data within section
        while True:
            line = f.readline()
            if ("$End%s" % title) not in line: newF.write(line) 
            else: break

        f.seek(0); newF.seek(0) # restore file positions
        return newF

    def _seekForHeader(self, title, f):
        """
        Iterate through a file until we end up at the section header
        for `title`. Function has obvious side-effects on `f`.
        """
        while True:
            line = f.readline()
            if len(line) == 0:
                raise EOFError("No `%s' header found!" % title)
                break
            elif (("$%s" % title) not in line): continue
            else: break 

    def _vertexCoordsAndMap(self):
        """
        Extract vertex coordinates and mapping information from
        nx.genfromtxt-friendly file, generated in `_isolateData`.

        Returns both the vertex coordinates and the mapping information.
        Mapping information is stored in a 1xn array where n is the
        largest vertexID obtained from the gmesh file. This mapping
        array is subsequently used to transform element information.
        """
        gen = nx.genfromtxt(fname=self.nodesFile, skiprows=1)
        self.nodesFile.close()

        vertexCoords = gen[:, 1:] # strip out column 0
        vertexIDs    = gen[:, :1].flatten().astype(int)
        
        # `vertexToIdx`: gmsh-vertex ID -> `vertexCoords` index
        vertexToIdx = nx.empty(vertexIDs.max() + 1)
        vertexToIdx[vertexIDs] = nx.arange(len(vertexIDs))

        # transpose for FiPy, truncate for dimension
        return vertexCoords.transpose()[:self.coordDimensions], vertexToIdx

    def _parseElements(self, vertexMap):
        """
        Parses $Elements portion to build `facesToVertices` and `cellsToFaces`
        arrays.
        """
        def _extractFaces(faceLen, facesPerCell, cell):
            """
            Given `cell`, a cell in terms of vertices, returns an array of
            `facesPerCell` faces of length `faceLen` in terms of vertices.
            """
            faces = []
            for i in range(facesPerCell):
                aFace = []
                for j in range(faceLen):
                    aVertex = (i + j) % len(cell) # we may wrap
                    aFace.append(int(cell[aVertex]))
                faces.append(aFace)
            return faces

        def makeExtractFacesFnc(faceLen, facesPerCell):
            """ Wrapper for `_extractFaces`. """
            return lambda c: _extractFaces(faceLen, facesPerCell, c)

        def formatForFiPy(arr): return arr.swapaxes(0,1)[::-1]

        cellsToVertIDs = []
        els            = self.elemsFile.readlines()

        # read in Elements data from gmsh
        for element in els[1:]: # skip number-of-elems line
            currLineInts = [int(x) for x in element.split()]
            elemType     = currLineInts[1]
            numTags      = currLineInts[2]

            if elemType in self.numFacesForShape.keys():
                # 3 columns precede the tags
                cellsToVertIDs.append(currLineInts[(3+numTags):])
            else:
                continue # shape not recognized

        self.elemsFile.close() # tempfile trashed

        # translate gmsh vertex IDs to vertexCoords indices
        cellsToVertices = vertexMap[nx.array(cellsToVertIDs, dtype=int)]

        # a few scalers
        # ASSUMPTION: all elements are of the same shape
        faceLength      = self.dimensions # number of vertices in a face
        numCells        = len(cellsToVertices)
        facesPerCell    = self.numFacesForShape[elemType]
        currNumFaces    = 0

        # a few data structures and a function
        cellsToFaces    = nx.empty((numCells, facesPerCell), dtype=int)
        facesDict       = {}
        uniqueFaces     = []
        facesFromCell   = makeExtractFacesFnc(faceLength, facesPerCell)

        # we now build, explicitly, `cellsToFaces` and `uniqueFaces`,
        # the latter will result in `facesToVertices`.
        for cellIdx in range(numCells):
            cell  = cellsToVertices[cellIdx]
            faces = facesFromCell(cell)

            for faceIdx in range(facesPerCell):
                currFace = faces[faceIdx]
                keyStr   = ' '.join([str(x) for x in sorted(currFace)])
                # NB: currFace is sorted for the key as to spot duplicates

                if facesDict.has_key(keyStr):
                    cellsToFaces[cellIdx][faceIdx] = facesDict[keyStr]
                else: # new face
                    facesDict[keyStr] = currNumFaces
                    cellsToFaces[cellIdx][faceIdx] = currNumFaces
                    uniqueFaces.append(currFace)
                    currNumFaces += 1

        facesToVertices = nx.array(uniqueFaces, dtype=int)

        return formatForFiPy(facesToVertices), formatForFiPy(cellsToFaces)

class GmshImporter2D(mesh2D.Mesh2D):
    """
    """
    def __init__(self, arg, coordDimensions=2):
        self.mshFile = MshFile(arg, dimensions=2, 
                               coordDimensions=coordDimensions)
        self.verts   = self.mshFile.vertexCoords
        self.faces   = self.mshFile.facesToV
        self.cells   = self.mshFile.cellsToF
        mesh2D.Mesh2D.__init__(self, vertexCoords=self.verts,
                                     faceVertexIDs=self.faces,
                                     cellFaceIDs=self.cells)
    def getCellVolumes(self):
        return abs(mesh2D.Mesh2D.getCellVolumes(self))

class GmshImporter2DIn3DSpace(GmshImporter2D):
    def __init__(self, arg):
        GmshImporter2D.__init__(self, arg, coordDimensions=3)

class GmshImporter3D(mesh.Mesh):
    def __init__(self, arg):
        self.mshFile = MshFile(arg, dimensions=3)
        self.verts   = self.mshFile.vertexCoords
        self.faces   = self.mshFile.facesToV
        self.cells   = self.mshFile.cellsToF
        mesh.Mesh.__init__(self, vertexCoords=self.verts,
                                 faceVertexIDs=self.faces,
                                 cellFaceIDs=self.cells)
    
