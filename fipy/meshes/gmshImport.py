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

__docformat__ = 'restructuredtext'

from fipy.tools import numerix as nx
from fipy.tools import parallel
from fipy.tools import serial
from fipy.tools.decorators import getsetDeprecated
from fipy.meshes.mesh import Mesh
from fipy.meshes.mesh2D import Mesh2D
import sys
import os
import tempfile

DEBUG = 0


def parprint(str):
    if DEBUG:
        if parallel.procID == 0:
            print >> sys.stderr, str

class GmshException(Exception):
    pass

class MshFile:
    """
    Class responsible for parsing a Gmsh file and then readying
    its contents for use by a `Mesh` constructor. 
    
    Can handle a partitioned mesh based on `parallel.Nproc`. If partitioning,
    the msh file must either be previously partitioned with the number of
    partitions matching `Nproc`, or the mesh must be specified with a geo file
    or multiline string.

    Does not support gmsh versions < 2. If partitioning, gmsh
    version must be >= 2.5.

    TODO: Refactor face extraction functions.
    """
    def __init__(self, filename, 
                       dimensions, 
                       coordDimensions=None,
                       communicator=parallel,
                       order=1):
        """
        Isolate relevant data into two files, store in 
        `self.nodesFile` for $Nodes,
        `self.elemsFile` for $Elements. 

        :Parameters:
          - `filename`: a string indicating gmsh output file
          - `dimensions`: an integer indicating dimension of mesh
          - `coordDimension`: an integer indicating dimension of shapes
        """
        
        self.communicator    = communicator
        self.coordDimensions = coordDimensions or dimensions
        self.dimensions      = dimensions

        if order > 1:
            self.communicator = serial

        # much special-casing based on gmsh version
        gmshVersion = self._gmshVersion
        if gmshVersion < 2.0:
            errStr = "Gmsh version must be >= 2.0."
            raise EnvironmentError(errStr)
        if self.communicator.Nproc > 1:
            if gmshVersion < 2.5:
                import warnings
                warnstr = "Cannot partition with Gmsh version < 2.5. " \
                           + "Reverting to serial."
                warnings.warn(warnstr, RuntimeWarning, stacklevel=2)
                self.communicator = serial

                gmshFlags = "-%d" % (self.dimensions)
            else: # gmsh version is adequate for partitioning
                gmshFlags = "-%d -part %d" % (self.dimensions,
                                              self.communicator.Nproc)
        else: # we're running serial
            gmshFlags = "-%d" % (self.dimensions)
        
        gmshFlags += " -format msh"

        self.filename, self.gmshOutput = self._parseFilename(filename, 
                                                             gmshFlags)

        # we need a conditional here so we don't pick up 2D shapes in 3D
        if dimensions == 2: 
            self.numFacesForShape = {2: 3, # triangle:   3 sides
                                     3: 4} # quadrangle: 4 sides
        else: # 3D
            self.numFacesForShape = {4: 4, # tet:        4 sides
                                     5: 6} # hexahedron: 6 sides

        self.faceLenForShape = {2: 2, # triangle:   2 verts per face
                                3: 2, # quadrangle: 2 verts per face
                                4: 3, # tet:        3 verts per face
                                5: 4} # hexahedron: 4 verts per face

        f = open(self.filename, "r") # open the msh file

        self.version, self.fileType, self.dataSize = self._getMetaData(f)
        self.nodesFile = self._isolateData("Nodes", f)
        self.elemsFile = self._isolateData("Elements", f)

    def _parseFilename(self, fname, gmshFlags):
        """
        If we're being passed a .msh file, leave it be. Otherwise,
        we've gotta compile a .msh file from either (i) a .geo file, 
        or (ii) a gmsh script passed as a string.
        """
        import subprocess as subp
        lowerFname = fname.lower()
        if '.msh' in lowerFname:
            return fname, None
        else:
            if '.geo' in lowerFname or '.gmsh' in lowerFname:
                geoFile = fname
            else: # fname must be a full script, not a file
                (f, geoFile) = tempfile.mkstemp('.geo')
                file = open(geoFile, 'w')
                file.writelines(fname)
                file.close(); os.close(f)

            (f, mshFile) = tempfile.mkstemp('.msh')
            gmshout = subp.Popen("gmsh %s %s -o %s" \
                      % (geoFile, gmshFlags, mshFile),
                      stdout=subp.PIPE, shell=True).stdout.readlines()
            parprint("gmsh out: %s" % gmshout)
            os.close(f)

            return mshFile, gmshout
         
    def _getMetaData(self, f):
        """
        Extracts gmshVersion, file-type, and data-size in that
        order.
        """
        self._seekForHeader("MeshFormat", f)
        metaData = f.readline().split()
        f.seek(0)
        return [float(x) for x in metaData]

    @property
    def _gmshVersion(self):
        """
        Enforce gmsh version to be either >= 2 or 2.5, based on Nproc.
        
        We can't trust the generated msh file for the correct version number, so
        we have to retrieve it from the gmsh binary.
        """
        import subprocess as subp
        import re

        verStr = "\n".join(subp.Popen("gmsh --version", 
            stderr=subp.PIPE, shell=True).stderr.readlines())

        m = re.search(r'\d+.\d+', verStr)

        if m:
            return float(m.group(0))
        else:
            return 0
     
    def _isolateData(self, title, f):
        """
        Gets all data between $[title] and $End[title], writes
        it out to its own file.
        """
        newF = tempfile.TemporaryFile()
        self._seekForHeader(title, f)
        
        # extract the actual data within section
        while True:
            line = f.readline()
            if ("$End%s" % title) not in line: 
                newF.write(line) 
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
            else: break # found header

    def _deriveCellsAndFaces(self, cellsToVertIDs, shapeTypes, numCells):
        """
        Uses element information obtained from `_parseElementFile` to deliver
        `facesToVertices` and `cellsToFaces`.
        """

        def formatForFiPy(arr): return arr.swapaxes(0,1)[::-1]

        allShapes  = nx.unique(shapeTypes).tolist()
        maxFaces   = max([self.numFacesForShape[x] for x in allShapes])

        # `cellsToFaces` must be padded with -1; see mesh.py
        currNumFaces = 0
        cellsToFaces = nx.ones((numCells, maxFaces)) * -1
        facesDict    = {}
        uniqueFaces  = []

        # we now build `cellsToFaces` and `uniqueFaces`,
        # the latter will result in `facesToVertices`.
        for cellIdx in range(numCells):
            shapeType    = shapeTypes[cellIdx]
            faceLength   = self.faceLenForShape[shapeType]
            cell         = cellsToVertIDs[cellIdx]
            facesPerCell = self.numFacesForShape[shapeType]

            if shapeType == 5: # we need to special case for hexahedron
                faces = self._extractHexahedronFaces(cell)
            else:
                faces = self._extractFaces(faceLength, facesPerCell, cell)

            for faceIdx in range(facesPerCell):
                # NB: currFace is sorted for the key to spot duplicates
                currFace = faces[faceIdx]
                keyStr   = ' '.join([str(x) for x in sorted(currFace)])

                if facesDict.has_key(keyStr):
                    cellsToFaces[cellIdx][faceIdx] = facesDict[keyStr]
                else: # new face
                    facesDict[keyStr] = currNumFaces
                    cellsToFaces[cellIdx][faceIdx] = currNumFaces
                    uniqueFaces.append(currFace)
                    currNumFaces += 1

        facesToVertices = nx.array(uniqueFaces, dtype=int)

        return formatForFiPy(facesToVertices), formatForFiPy(cellsToFaces)
 
    def _translateVertIDToIdx(self, cellsToVertIDs, vertexMap):
        """
        Translates cellToIds from Gmsh output IDs to `vertexCoords`
        indices. 
        
        NB: Takes in Python array, outputs numpy array.
        """
        cellsToVertIdxs = []

        # translate gmsh vertex IDs to vertexCoords indices
        for cell in cellsToVertIDs:
            vertIndices = vertexMap[nx.array(cell)]
            cellsToVertIdxs.append(vertIndices)

        return cellsToVertIdxs
     
    def _extractFaces(self, faceLen, facesPerCell, cell):
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

    def _extractHexahedronFaces(self, cell):
        """
        SPECIAL CASE: return faces for a hexahedron cell.
        """
        def orderingToFace(vertList):
            aFace = []
            for i in vertList:
                aFace.append(int(cell[i]))
            return aFace

        # six orderings for six faces
        faces = []
        orderings = [[0, 1, 2, 3], # ordering of vertices gleaned from
                     [4, 5, 6, 7], # a one-cube Grid3D example
                     [0, 1, 5, 4],
                     [3, 2, 6, 7],
                     [0, 3, 7, 4],
                     [1, 2, 6, 5]]

        for o in orderings:
            faces.append(orderingToFace(o))

        return faces

    def buildMeshData(self):
        """
        0. Build cellsToVertices
        1. Recover needed vertexCoords and mapping from file using 
           cellsToVertices
        2. Build cellsToVertIDs proper from vertexCoords and vertex map
        3. Build faces
        4. Build cellsToFaces

        Returns vertexCoords, facesToVertexID, cellsToFaceID, 
                cellGlobalIDMap, ghostCellGlobalIDMap.
        """
     
        parprint("Parsing elements.")
        cellDataDict, ghostDataDict = self._parseElementFile()
        
        cellsToGmshVerts = cellDataDict['cells']  + ghostDataDict['cells']
        numCellsTotal    = cellDataDict['num']    + ghostDataDict['num']
        allShapeTypes    = cellDataDict['shapes'] + ghostDataDict['shapes']
        allShapeTypes    = nx.array(allShapeTypes)
        allShapeTypes    = nx.delete(allShapeTypes, nx.s_[numCellsTotal:])

        if numCellsTotal < 1:
            errStr = "Gmsh hasn't produced any cells! Check your Gmsh code."
            errStr += "\n\nGmsh output:\n%s" % "".join(self.gmshOutput).rstrip()
            raise GmshException(errStr)

        parprint("Recovering coords.")
        parprint("numcells %d" % numCellsTotal)
        vertexCoords, vertIDtoIdx = self._vertexCoordsAndMap(cellsToGmshVerts)

        # translate Gmsh IDs to `vertexCoord` indices
        cellsToVertIDs = self._translateVertIDToIdx(cellsToGmshVerts,
                                                    vertIDtoIdx)

        parprint("Building cells and faces.")
        facesToV, cellsToF = self._deriveCellsAndFaces(cellsToVertIDs,
                                                       allShapeTypes,
                                                       numCellsTotal)
        parprint("Done with cells and faces.")
        return vertexCoords, facesToV, cellsToF, \
               cellDataDict['idmap'], ghostDataDict['idmap']

    def _vertexCoordsAndMap(self, cellsToGmshVerts):
        """
        Returns `vertexCoords` and mapping from Gmsh ID to `vertexCoords`
        indices (same as in MshFile). 
        
        Unlike parent, doesn't use genfromtxt
        because we want to avoid loading the entire msh file into memory.
        """
        allVerts     = [v for c in cellsToGmshVerts for v in c] # flatten
        allVerts     = nx.unique(nx.array(allVerts, dtype=int)) # remove dups
        allVerts     = nx.sort(allVerts)
        maxVertIdx   = allVerts[-1] + 1 # add one to offset zero
        vertGIDtoIdx = nx.ones(maxVertIdx) * -1 # gmsh ID -> vertexCoords idx
        vertexCoords = nx.empty((len(allVerts), self.coordDimensions))
        nodeCount    = 0

        # establish map. This works because allVerts is a sorted set.
        vertGIDtoIdx[allVerts] = nx.arange(len(allVerts))

        self.nodesFile.readline() # skip number of nodes

        # now we walk through node file with a sorted unique list of vertices
        # in hand. When we encounter 0th element in `allVerts`, save it 
        # to `vertexCoords` then pop its ID off `allVerts`.
        node = self.nodesFile.readline()
        while node != "":
            line   = node.split()
            nodeID = int(line[0])

            if nodeID == allVerts[nodeCount]:
                newVert = [float(x) for x in line[1:self.coordDimensions+1]]
                vertexCoords[nodeCount,:] = nx.array(newVert)
                nodeCount += 1

            if len(allVerts) == nodeCount: 
                break

            node = self.nodesFile.readline()

        # transpose for FiPy
        transCoords = vertexCoords.swapaxes(0,1)
        return transCoords, vertGIDtoIdx

    def _parseElementFile(self):
        """
        Return two dicts, the first for non-ghost cells and the second for
        ghost cells. Each dict contains 
            "cells": A cellToVertID Python array,
            "shapes": A shapeTypes Python array,
            "num": A scalar, number of cells
            "idmap": A Python array which maps vertexCoords idx -> global ID

        All nastiness concerning ghost cell
        calculation is consolidated here: if we were ever to need to CALCULATE
        GHOST CELLS OURSELVES, the only code we'd have to change is in here.
        """
        def _addCell(cellDict, currLine, elType, IDoffset):
            """
            Bookkeeping for cells. Declared as own function for generality.
            Curried later in ghost/non-ghost specific lambdas.
            Just full of side-effects. Sorry, Abelson/Sussman.
            """
            cellDict['cells'].append(currLine[(numTags+3):])
            cellDict['shapes'].append(elType)
            cellDict['idmap'].append(currLine[0] - IDoffset)
            cellDict['num'] += 1

        seenParts  = [] # to make sure Nproc/mshfile are using same num. proc.
        cellsData  = {'cells':  [],
                      'shapes': [],
                      'num':     0,
                      'idmap':  []} # vertexCoords idx -> gmsh ID (global ID)

        ghostsData = {'cells':  [],
                      'shapes': [],
                      'num':     0,
                      'idmap':  []} # vertexCoords idx -> gmsh ID (global ID)

        IDOffset        = -1 # this will be subtractd from gmsh ID to obtain
                             # global ID
        pid             = self.communicator.procID + 1

        addCell      = lambda c,e: _addCell(cellsData, c, e, IDOffset)
        addGhostCell = lambda c,e: _addCell(ghostsData, c, e, IDOffset)

        self.elemsFile.readline() # skip number of elements
        for el in self.elemsFile:
            currLineInts = [int(x) for x in el.split()]
            elemType     = currLineInts[1]

            if elemType in self.numFacesForShape.keys():
                if IDOffset == -1: # if first valid shape
                    IDOffset = currLineInts[0]

                numTags = currLineInts[2]
                tags    = currLineInts[3:(3+numTags)]
                # if we're collecting ghost cells
                if self.communicator.Nproc > 1: 
                    tags.reverse() # any ghost cell IDs come first, then part ID

                    if tags[0] < 0: # if this is someone's ghost cell
                        while tags[0] < 0: # while we have ghost IDs
                            if (tags[0] * -1) == pid: 
                                addGhostCell(currLineInts, elemType)
                                break   
                            tags.pop(0) # try the next ID
                # el is in this processor's partition OR we collect all cells
                if (self.communicator.Nproc == 1) or (tags[0] == pid):
                    addCell(currLineInts, elemType)
                
        self.elemsFile.close() # tempfile trashed

        return cellsData, ghostsData

class Gmsh2D(Mesh2D):
    def __init__(self, 
                 arg, 
                 coordDimensions=2, 
                 communicator=parallel, 
                 order=1):
        self.communicator = communicator 
        self.mshFile  = MshFile(arg, 
                                dimensions=2, 
                                coordDimensions=coordDimensions,
                                communicator=communicator,
                                order=order)
        (verts,
        faces,
        cells,
        self.cellGlobalIDs, 
        self.gCellGlobalIDs) = self.mshFile.buildMeshData()

        if self.communicator.Nproc > 1:
            self.globalNumberOfCells = self.communicator.sumAll(len(self.cellGlobalIDs))
            parprint("  I'm solving with %d cells total." % self.globalNumberOfCells)
            parprint("  Got global number of cells")

        Mesh2D.__init__(self, vertexCoords=verts,
                              faceVertexIDs=faces,
                              cellFaceIDs=cells)
        parprint("Exiting Gmsh2D")

    def __setstate__(self, dict):
        Mesh2D.__init__(self, **dict)
        self.cellGlobalIDs = list(nx.arange(self.cellFaceIDs.shape[-1]))
        self.gCellGlobalIDs = []
        self.communicator = serial
        self.mshFile = None
    
    @getsetDeprecated
    def _getGlobalNonOverlappingCellIDs(self):
        return self._globalNonOverlappingCellIDs

    @property
    def _globalNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Does not include the IDs of boundary cells.

        E.g., would return [0, 1, 4, 5] for mesh A

            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """ 
        return nx.array(self.cellGlobalIDs)

    @getsetDeprecated
    def _getGlobalOverlappingCellIDs(self):
        return self._globalOverlappingCellIDs

    @property
    def _globalOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 4, 5, 6] for mesh A

            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """ 
        return nx.array(self.cellGlobalIDs + self.gCellGlobalIDs)

    @getsetDeprecated
    def _getLocalNonOverlappingCellIDs(self):
        return self._localNonOverlappingCellIDs

    @property
    def _localNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

            A        B
        ------------------
        | 3 | 4 || 4 | 5 |
        ------------------
        | 0 | 1 || 1 | 2 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return nx.arange(len(self.cellGlobalIDs))

    @getsetDeprecated
    def _getLocalOverlappingCellIDs(self):
        return self._localOverlappingCellIDs

    @property
    def _localOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3, 4, 5] for mesh A

            A        B
        ------------------
        | 3 | 4 || 5 |   |
        ------------------
        | 0 | 1 || 2 |   |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return nx.arange(len(self.cellGlobalIDs) 
                         + len(self.gCellGlobalIDs))
    
    def _test(self):
        """
        First, we'll test GmshImporter2D on a small circle with triangular
        cells.

        >>> circ = Gmsh2D('''
        ... cellSize = 1; 
        ... radius   = 0.25; 
        ... Point(1) = {0, 0, 0, cellSize}; 
        ... Point(2) = {-radius, 0, 0, cellSize}; 
        ... Point(3) = {0, radius, 0, cellSize}; 
        ... Point(4) = {radius, 0, 0, cellSize}; 
        ... Point(5) = {0, -radius, 0, cellSize}; 
        ... Circle(6) = {2, 1, 3}; 
        ... Circle(7) = {3, 1, 4}; 
        ... Circle(8) = {4, 1, 5}; 
        ... Circle(9) = {5, 1, 2}; 
        ... Line Loop(10) = {6, 7, 8, 9}; 
        ... Plane Surface(11) = {10}; 
        ... ''')

        >>> print circ.cellVolumes[0] > 0
        True

        Now we'll test GmshImporter2D again, but on a rectangle.

        >>> rect = Gmsh2D('''
        ... cellSize = 0.5;
        ... radius   = 10;
        ... Point(2) = {-radius, radius, 0, cellSize};
        ... Point(3) = {radius, radius, 0, cellSize};
        ... Point(4) = {radius, -radius, 0, cellSize};
        ... Point(5) = {-radius, -radius, 0, cellSize};
        ... Line(6) = {2, 3};
        ... Line(7) = {3, 4};
        ... Line(8) = {4, 5};
        ... Line(9) = {5, 2};
        ... Line Loop(10) = {6, 7, 8, 9};
        ... Plane Surface(11) = {10};
        ... ''')

        >>> print rect.cellVolumes[0] > 0
        True

        Testing multiple shape types within a mesh;

        >>> circle = Gmsh2D('''
        ... cellSize = 0.05;
        ... radius = 1;
        ... Point(1) = {0, 0, 0, cellSize};
        ... Point(2) = {-radius, 0, 0, cellSize};
        ... Point(3) = {0, radius, 0, cellSize};
        ... Point(4) = {radius, 0, 0, cellSize};
        ... Point(5) = {0, -radius, 0, cellSize};
        ... Circle(6) = {2, 1, 3};
        ... Circle(7) = {3, 1, 4};
        ... Circle(8) = {4, 1, 5};
        ... Circle(9) = {5, 1, 2};
        ... Line Loop(10) = {6, 7, 8, 9};
        ... Plane Surface(11) = {10};
        ... Recombine Surface{11};
        ... ''')

        >>> print circle.cellVolumes[0] > 0
        True
        
        >>> from fipy.tools import dump
        >>> from fipy import parallel
        >>> f, tempfile = dump.write(circle)
        >>> pickle_circle = dump.read(tempfile, f)

        >>> print parallel.Nproc > 1 or (pickle_circle.cellVolumes == circle.cellVolumes).all()
        True
        
        >>> print parallel.Nproc > 1 or (pickle_circle._globalOverlappingCellIDs == circle._globalOverlappingCellIDs).all()
        True

        >>> cmd = "Point(1) = {0, 0, 0, 0.05};"

        >>> Gmsh2D(cmd) #doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        GmshException: Gmsh hasn't produced any cells! Check your Gmsh code.
        """

class Gmsh2DIn3DSpace(Gmsh2D):
    def __init__(self, arg, communicator=parallel, order=1):
        Gmsh2D.__init__(self, 
                        arg, 
                        coordDimensions=3, 
                        communicator=communicator,
                        order=order)

    def _test(self):
        """
        Stolen from the cahnHilliard sphere example.

        >>> sphere = Gmsh2DIn3DSpace('''
        ... radius = 5.0;
        ... cellSize = 0.3;
        ...
        ... // create inner 1/8 shell
        ... Point(1) = {0, 0, 0, cellSize};
        ... Point(2) = {-radius, 0, 0, cellSize};
        ... Point(3) = {0, radius, 0, cellSize};
        ... Point(4) = {0, 0, radius, cellSize};
        ... Circle(1) = {2, 1, 3};
        ... Circle(2) = {4, 1, 2};
        ... Circle(3) = {4, 1, 3};
        ... Line Loop(1) = {1, -3, 2} ;
        ... Ruled Surface(1) = {1};
        ...
        ... // create remaining 7/8 inner shells
        ... t1[] = Rotate {{0,0,1},{0,0,0},Pi/2}
        ... {Duplicata{Surface{1};}};
        ... t2[] = Rotate {{0,0,1},{0,0,0},Pi}
        ... {Duplicata{Surface{1};}};
        ... t3[] = Rotate {{0,0,1},{0,0,0},Pi*3/2}
        ... {Duplicata{Surface{1};}};
        ... t4[] = Rotate {{0,1,0},{0,0,0},-Pi/2}
        ... {Duplicata{Surface{1};}};
        ... t5[] = Rotate {{0,0,1},{0,0,0},Pi/2}
        ... {Duplicata{Surface{t4[0]};}};
        ... t6[] = Rotate {{0,0,1},{0,0,0},Pi}
        ... {Duplicata{Surface{t4[0]};}};
        ... t7[] = Rotate {{0,0,1},{0,0,0},Pi*3/2}
        ... {Duplicata{Surface{t4[0]};}};
        ...
        ... // create entire inner and outer shell
        ... Surface
        ... Loop(100)={1,t1[0],t2[0],t3[0],t7[0],t4[0],t5[0],t6[0]};
        ... ''').extrude(extrudeFunc=lambda r: 1.1 * r)

        >>> print sphere.cellVolumes[0] > 0
        True

        >>> from fipy.tools import dump
        >>> from fipy import parallel
        >>> f, tempfile = dump.write(sphere)
        >>> pickle_sphere = dump.read(tempfile, f)

        >>> print parallel.Nproc > 1 or (pickle_sphere.cellVolumes == sphere.cellVolumes).all()
        True
        
        >>> print parallel.Nproc > 1 or (pickle_sphere._globalOverlappingCellIDs == sphere._globalOverlappingCellIDs).all()
        True
        """

class Gmsh3D(Mesh):
    def __init__(self, arg, communicator=parallel, order=1):
        self.mshFile  = MshFile(arg, 
                                dimensions=3, 
                                communicator=communicator,
                                order=order)

        (verts,
        faces,
        cells,
        self.cellGlobalIDs,
        self.gCellGlobalIDs) = self.mshFile.buildMeshData()

        Mesh.__init__(self, vertexCoords=verts,
                            faceVertexIDs=faces,
                            cellFaceIDs=cells)

        self.communicator = communicator

        if self.communicator.Nproc > 1:
            self.globalNumberOfCells = self.communicator.sumAll(len(self.cellGlobalIDs))

    def __setstate__(self, dict):
        Mesh.__init__(self, **dict)
        self.cellGlobalIDs = list(nx.arange(self.cellFaceIDs.shape[-1]))
        self.gCellGlobalIDs = []
        self.communicator = serial
        self.mshFile = None

    @getsetDeprecated
    def _getGlobalNonOverlappingCellIDs(self):
        return self._globalNonOverlappingCellIDs

    @property
    def _globalNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Does not include the IDs of boundary cells.

        E.g., would return [0, 1, 4, 5] for mesh A

            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """ 
        return nx.array(self.cellGlobalIDs)

    @getsetDeprecated
    def _getGlobalOverlappingCellIDs(self):
        return self._globalOverlappingCellIDs

    @property
    def _globalOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 4, 5, 6] for mesh A

            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """ 
        return nx.array(self.cellGlobalIDs + self.gCellGlobalIDs)

    @getsetDeprecated
    def _getLocalNonOverlappingCellIDs(self):
        return self._localNonOverlappingCellIDs

    @property
    def _localNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

            A        B
        ------------------
        | 3 | 4 || 4 | 5 |
        ------------------
        | 0 | 1 || 1 | 2 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return nx.arange(len(self.cellGlobalIDs))

    @getsetDeprecated
    def _getLocalOverlappingCellIDs(self):
        return self._localOverlappingCellIDs

    @property
    def _localOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3, 4, 5] for mesh A

            A        B
        ------------------
        | 3 | 4 || 5 |   |
        ------------------
        | 0 | 1 || 2 |   |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return nx.arange(len(self.cellGlobalIDs) 
                         + len(self.gCellGlobalIDs))
     
    def _test(self):
        """
        >>> prism = Gmsh3D('''
        ... cellSize = 0.5;
        ... Len = 2;
        ... Hei = 1;
        ... Wid = 1;
        ...
        ... Point(1) = {0, 0, 0, cellSize};
        ... Point(2) = {0, 0, Wid, cellSize};
        ... Point(3) = {0, Hei, Wid, cellSize};
        ... Point(4) = {0, Hei, 0, cellSize};
        ...
        ... Point(5) = {Len, 0, 0, cellSize};
        ... Point(6) = {Len, 0, Wid, cellSize};
        ... Point(7) = {Len, Hei, Wid, cellSize};
        ... Point(8) = {Len, Hei, 0, cellSize};
        ...
        ... Line(9)  = {1, 2};
        ... Line(10) = {2, 3};
        ... Line(11) = {3, 4};
        ... Line(12) = {4, 1};
        ...
        ... Line(13) = {5, 6};
        ... Line(14) = {6, 7};
        ... Line(15) = {7, 8};
        ... Line(16) = {8, 5};
        ...
        ... Line(17) = {1, 5};
        ... Line(18) = {2, 6};
        ... Line(19) = {3, 7};
        ... Line(20) = {4, 8};
        ...
        ... Line Loop(21) = {9, 10, 11, 12};
        ... Line Loop(22) = {13, 14, 15, 16};
        ... Line Loop(23) = {17, -16, -20, 12};
        ... Line Loop(24) = {13, -18, -9, 17};
        ... Line Loop(25) = {18, 14, -19, -10};
        ... Line Loop(26) = {-19, 11, 20, -15};
        ...
        ... Plane Surface(27) = {21};
        ... Plane Surface(28) = {22};
        ... Plane Surface(29) = {23};
        ... Plane Surface(30) = {24};
        ... Plane Surface(31) = {25};
        ... Plane Surface(32) = {26};
        ...
        ... Surface Loop(33) = {27, 28, 29, 30, 31, 32};
        ...
        ... Volume(34) = {33};
        ... ''')

        >>> print prism.cellVolumes[0] > 0
        True
        
        >>> from fipy.tools import dump
        >>> from fipy import parallel
        >>> f, tempfile = dump.write(prism)
        >>> pickle_prism = dump.read(tempfile, f)

        >>> print parallel.Nproc > 1 or (pickle_prism.cellVolumes == prism.cellVolumes).all()
        True
        
        >>> print parallel.Nproc > 1 or (pickle_prism._globalOverlappingCellIDs == prism._globalOverlappingCellIDs).all()
        True
        """

class GmshGrid2D(Gmsh2D):
    """Should serve as a drop-in replacement for Grid2D."""
    def __init__(self, dx=1., dy=1., nx=1, ny=None, 
                 coordDimensions=2, communicator=parallel, order=1):
        self.dx = dx
        self.dy = dy or dx
        self.nx = nx
        self.ny = ny or nx

        arg = self._makeGridGeo(self.dx, self.dy, self.nx, self.ny)

        Gmsh2D.__init__(self, arg, coordDimensions, communicator, order)
    
    @getsetDeprecated
    def _getMeshSpacing(self):
        return self._meshSpacing

    @property
    def _meshSpacing(self):
        return nx.array((self.dx,self.dy))[...,nx.newaxis]

    def _makeGridGeo(self, dx, dy, nx, ny):
        height = ny * dy
        width  = nx * dx
        numLayers = int(ny / dy)

        # kludge: must offset cellSize by `eps` to work properly
        eps = float(dx)/(nx * 10) 

        return """
            ny       = %(ny)g;
            nx       = %(nx)g;
            cellSize = %(dx)g - %(eps)g;
            height   = %(height)g;
            width    = %(width)g;

            Point(1) = {0, 0, 0, cellSize};
            Point(2) = {width, 0, 0, cellSize};
            Line(3) = {1, 2};
            Extrude{0, height, 0} {
                    Line{3}; Layers{ ny }; Recombine;
            }
            """ % locals()   

    def _test(self):
        """
        Here we do some rudimentary comparisons between GmshGrid and Grid, just
        to ensure they're performing similarly.

        >>> from fipy import *

        >>> yogmsh = GmshGrid2D(dx=5, dy=5, nx=5, ny=5, communicator=serial)

        >>> yogrid = Grid2D(dx=5, dy=5, nx=5, ny=5, communicator=serial)

        >>> numerix.allclose(yogmsh._faceAreas, yogrid._faceAreas)
        True

        >>> yogmsh.cellCenters.value.size == yogrid.cellCenters.value.size
        True

        >>> mesh = GmshGrid2D(nx=2, ny=2)

        >>> mesh.numberOfCells == 4
        True

        >>> len(mesh.faceCenters[0]) == 12
        True
        """


class GmshGrid3D(Gmsh3D):
    """Should serve as a drop-in replacement for Grid3D."""
    def __init__(self, dx=1., dy=1., dz=1., nx=1, ny=None, nz=None,
                 communicator=parallel, order=1):
        self.dx = dx
        self.dy = dy or dx
        self.dz = dz or dx

        self.nx = nx
        self.ny = ny or nx
        self.nz = nz or nx

        arg = self._makeGridGeo(self.dx, self.dy, self.dz,
                                self.nx, self.ny, self.nz)

        Gmsh3D.__init__(self, arg, communicator=communicator, order=order)
    
    @getsetDeprecated
    def _getMeshSpacing(self):
        return self._meshSpacing

    @property
    def _meshSpacing(self):
        return nx.array((self.dx,self.dy,self.dz))[...,nx.newaxis]
 
    def _makeGridGeo(self, dx, dy, dz, nx, ny, nz):
        height = ny * dy
        width  = nx * dx
        depth  = nz * dz

        eps = float(dx)/(nx * 10)
        
        return """
            ny       = %(ny)g;
            nx       = %(nx)g;
            nz       = %(nz)g;
            cellSize = %(dx)g - %(eps)g;
            height   = %(height)g;
            width    = %(width)g;
            depth    = %(depth)g;

            Point(1) = {0, 0, 0, cellSize};
            Point(2) = {width, 0, 0, cellSize};
            Line(3) = {1, 2};
            out[] = Extrude{0, height, 0} {
                Line{3}; Layers{ ny }; Recombine;
            };
            Extrude{0, 0, depth} {
                Surface{ out[1] }; Layers{ nz }; Recombine;
            }
            """ % locals()   

    def _test(self):
        """
        Here we do some rudimentary comparisons between GmshGrid and Grid, just
        to ensure they're performing similarly.

        >>> from fipy import *
        >>> from fipy.tools import numerix as nx

        >>> yogmsh = GmshGrid3D(dx=5, dy=5, dz=5, nx=5, ny=5, nz=5,
        ...                     communicator=serial)

        >>> yogrid = Grid3D(dx=5, dy=5, dz=5, nx=5, ny=5, nz=5,
        ...                 communicator=serial)

        >>> yogmsh.cellCenters.value.size == yogrid.cellCenters.value.size
        True

        >>> numerix.allclose(yogmsh._faceAreas, yogrid._faceAreas)
        True

        >>> numerix.allclose(yogmsh._faceAreas, yogrid._faceAreas)
        True

        >>> mesh = GmshGrid3D(nx=2, ny=2, nz=2)

        >>> ccs = [[ 0.5,  0.5,  0.5,  0.5,  1.5,  1.5,  1.5,  1.5],
        ...    [ 0.5,  0.5,  1.5,  1.5,  0.5,  0.5,  1.5,  1.5],
        ...    [ 0.5,  1.5,  0.5,  1.5,  0.5,  1.5,  0.5,  1.5]]

        >>> len(mesh.cellCenters.value[0]) == 8
        True

        >>> faceAreas = [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
        ...           1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
        ...           1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.]

        >>> nx.allclose(mesh._faceAreas, faceAreas)
        True

        >>> cellAreas = [[ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
        ...            [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
        ...            [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
        ...            [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
        ...            [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
        ...            [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.]]

        >>> nx.allclose(mesh._cellAreas, cellAreas)
        True
        """
 
def deprecation(old, new):
    import warnings
    warnings.warn("%s has been replaced by %s." % (old, new), 
                  DeprecationWarning, stacklevel=3)

class GmshImporter2D(Gmsh2D):
    def __init__(self, arg, coordDimensions=2):
        deprecation("GmshImporter2D", "Gmsh2D")
        Gmsh2D.__init__(self, arg, coordDimensions=coordDimensions)

class GmshImporter2DIn3DSpace(Gmsh2DIn3DSpace):
    def __init__(self, arg):
        deprecation("GmshImporter2DIn3DSpace", "Gmsh2DIn3DSpace")
        Gmsh2DIn3DSpace.__init__(self, arg)

class GmshImporter3D(Gmsh3D):
    def __init__(self, arg):
        deprecation("GmshImporter3D", "Gmsh3D")
        Gmsh3D.__init__(self, arg)
    
def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
