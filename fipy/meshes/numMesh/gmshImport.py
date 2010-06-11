r"""
Does not support gmsh versions <= 2.
"""
from fipy.tools import numerix as np
import mesh
import mesh2D
import os

class MshFile:
    """
    Class responsible for parsing a Gmsh file and then readying
    its contents for use by a `Mesh` constructor.
    """
    def __init__(self, filename, dimension, coordDimensions=None):
        """
        Isolates relevant data into two files, stores in 
        `self.nodesFilename` for $Nodes,
        `self.elemsFilename` for $Elements.

        :Parameters:
          - `filename`: a string indicating gmsh output file
          - `dimension`: an integer indicating dimension of mesh
          - `coordDimension`: an integer indicating dimension of shapes

        TODO: Use tempfiles.
        """
        
        self.coordDimensions  = coordDimensions or dimension
        self.dimension        = dimension
        self.filename         = self._parseFilename(filename)
        self.numFacesForShape = {2: 3, # triangle:   3 sides
                                 3: 4, # quadrangle: 4 sides
                                 4: 4} # tet:        4 sides

        f = open(self.filename, "r") # open the msh file

        # five lines of muck here allow most of the other methods to be
        # free of side-effects.
        self.version, self.fileType, self.dataSize = self._getMetaData(f)
        self.nodesFilename = self._isolateData("Nodes", f)
        self.elemsFilename = self._isolateData("Elements", f)
        self.vertexCoords, self.vertexMap = self._vertexCoordsAndMap()
        self.facesToV, self.cellsToF  = self._parseElements(self.vertexMap)

    def _parseFilename(self, fname):
        """
        If we're being passed a .msh file, leave it be. Otherwise,
        we've gotta compile a .msh file from either (i) a .geo file, 
        or (ii) a gmsh script passed as a string.
        """
        lowerFname = fname.lower()
        if '.msh' in lowerFname:
            return fname
        else:
            import tempfile
            if '.geo' in lowerFname or '.gmsh' in lowerFname:
                geoFile = fname
            else: # fname must be a full script, not a file
                (f, geoFile) = tempfile.mkstemp('.geo')
                file = open(geoFile, 'w')
                file.writelines(fname)
                file.close(); os.close(f)

            (f, mshFile) = tempfile.mkstemp('.msh')
            os.system('gmsh %s -2 -v 0 -format msh -o %s' \
                      % (geoFile, mshFile))
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
        newFilename = "%s.%s" % (self.filename, title)
        newF = open(newFilename, "w")

        # seek to section header
        self._seekForHeader(title, f)
        
        # extract the actual data within section
        while True:
            line = f.readline()
            if ("$End%s" % title) not in line: newF.write(line) 
            else: break

        f.seek(0); newF.close() # restore file position, close up
        return newFilename

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
        np.genfromtxt-friendly file, generated in `_isolateData`.

        Returns both the vertex coordinates and the mapping information.
        Mapping information is stored in a 1xn array where n is the
        largest vertexID obtained from the gmesh file. This mapping
        array is subsequently used to transform element information.
        """
        gen = np.genfromtxt(fname=self.nodesFilename, skiprows=1)
        vertexCoords = gen[:, 1:] # strip out column 0
        vertexIDs    = gen[:, :1].flatten().astype(int)
        
        # `vertexToIdx`: gmsh-vertex ID -> `vertexCoords` index
        vertexToIdx = np.empty(vertexIDs.max() + 1)
        vertexToIdx[vertexIDs] = np.arange(len(vertexIDs))

        # transpose for FiPy, truncate for dimension
        return vertexCoords.transpose()[:self.dimension], vertexToIdx

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
        elsFile = open(self.elemsFilename, "r")
        els     = elsFile.readlines()

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

        elsFile.close()

        # translate gmsh vertex IDs to vertexCoords indices
        cellsToVertices = vertexMap[np.array(cellsToVertIDs, dtype=int)]

        # a few scalers 
        faceLength      = self.dimension # number of vertices in a face
        numCells        = len(cellsToVertices)
        facesPerCell    = self.numFacesForShape[elemType]
        currNumFaces    = 0

        # a few data structures and a function
        cellsToFaces    = np.empty((numCells, facesPerCell), dtype=int)
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

        facesToVertices = np.array(uniqueFaces, dtype=int)

        return formatForFiPy(facesToVertices), formatForFiPy(cellsToFaces)

class GmshImporter2D(mesh2D.Mesh2D):
    """
    """
    def __init__(self, arg, coordDimensions=2):
        mshFile = MshFile(arg, dimension=2)
        verts   = mshFile.vertexCoords
        faces   = mshFile.facesToV
        cells   = mshFile.cellsToF
        mesh2D.Mesh2D.__init__(self, vertexCoords=verts,
                                     faceVertexIDs=faces,
                                     cellFaceIDs=cells)
    def getCellVolumes(self):
        return abs(mesh.Mesh.getCellVolumes(self))


class GmshImporter2DIn3DSpace(GmshImporter2D):
    def __init__(self, arg):
        GmshImporter2D.__init__(self, arg, coordDimensions=3)

class GmshImporter3D(mesh.Mesh):
    def __init__(self, arg):
        mshFile = MshFile(arg, dimension=3)
        verts   = mshFile.vertexCoords
        faces   = mshFile.facesToV
        cells   = mshFile.cellsToF
        mesh.Mesh.__init__(self, vertexCoords=verts,
                                 faceVertexIDs=faces,
                                 cellFaceIDs=cells)
    
