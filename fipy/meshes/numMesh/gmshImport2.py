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
    its contents for import via numpy.genfromtxt.
    """
    def __init__(self, filename, dimension, coordDimensions=None):
        """
        Isolates relevant data into two files, stores in 
        `self.nodesFilename` for $Nodes,
        `self.elemsFilename` for $Elements.

        :Parameters:
          - `filename`: a string indicating gmsh output file

        TODO: Use tempfiles.
        """
        
        self.coordDimensions = coordDimensions or dimension
        self.dimension       = dimension
        self.filename        = self._parseFilename(filename)
        # 2 -> triangle, 3 -> quadrangle, 4 -> tetrahedron
        self.acceptedShapes  = (2, 3, 4)

        f = open(self.filename, "r") # open the msh file

        # five lines of muck here allow the other methods to be
        # free of side-effects
        self.version, self.fileType, self.dataSize = self._getMetaData(f)
        self.nodesFilename = self._isolateData("Nodes", f)
        self.elemsFilename = self._isolateData("Elements", f)
        self.vertexCoords, self.vertexMap = self._vertexCoordsAndMap()
        self.verticesOfCells = self._parseElements(self.vertexMap)
        self.faces, self.cells = self._deriveFacesAndCells(
                                   self.verticesOfCells
                                 )

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
        while "$MeshFormat" not in f.readline(): continue
        return [int(x) for x in f.readline().split()]

    def _isolateData(self, title, f):
        """
        Gets all data between $[title] and $End[title], writes
        it out to its own file.

        Very, very duct-tapey since it doesn't attempt to detect
        the beginning of a section; it will assume `title` is 
        the first section in `lines`.
        """
        newFilename = "%s.%s" % (self.filename, title)
        newF = open(newFilename, "w")

        # seek to section header
        while ("$%s" % title) not in f.readline(): continue

        # extract the actual data within section
        while True:
            line = f.readline()
            if ("$End%s" % title) not in line: newF.write(line) 
            else: break

        f.seek(0); newF.close() # restore file position, close up
        return newFilename

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
        Extracts and returns which nodes makeup which cells in
        the following format (example of triangular cells):

            cell0             celln
              -------      ------
                    |      |
        vertexIdx [ 0 ... 101 ]
        vertexIdx | 1 ... 102 |
        vertexIdx [ 2 ... 0   ]

        .. note:: A ``vertexIdx'' is an index for `vertexCoords`.
        """
        cellsToVertIDs = [] # will generate matrix described above
        elsFile = open(self.elemsFilename, "r")
        els     = elsFile.readlines()

        for element in els[1:]: # skip number-of-elems line
            currLineInts = [int(x) for x in element.split()]
            elemType     = currLineInts[1]
            numTags      = currLineInts[2]

            if elemType in self.acceptedShapes:
                cellsToVertIDs.append(
                    currLineInts[(3+numTags):] # 3 col. precede the tags
                )
            else:
                continue # shape not recognized

        elsFile.close()

        # transpose and then map vertex IDs to `vertexCoords` indices
        return vertexMap[np.array(cellsToVertIDs).transpose()]

    def _deriveFacesAndCells(self, verticesOfCells):
        """
        Returns two numpy arrays:
            `noDups`: faces in terms of vertices,
            `cellsToFaces`: cells in terms of faces.

        Unfortunately, it is most convenient to do this all 
        within one function.
        """
        numCells             = verticesOfCells.shape[-1]
        numVerticesInElement = len(verticesOfCells)
        cellFaceVertexIDs    = np.ones((self.dimension,
                                        numVerticesInElement,
                                        numCells)) * -1

        # build permutation matrix of valid faces
        for i in range(numVerticesInElement):
            cellFaceVertexIDs[:,i,:] = np.delete(verticesOfCells, i, 0)

        # reshape, sort so we can kill duplicates
        # NB: `sortedFaces` is the in the same ``order'' as
        # `cellFaceVertexIDs`, some of the pairs are just flipped
        cellFaceVertexIDs = cellFaceVertexIDs.reshape((self.dimension, -1))
        sortedFaces       = np.sort(cellFaceVertexIDs, axis=0)

        # filter duplicates with a hashmap while building `facesIdxMap`,
        # a mapping from indices into `cellFaceVertexIDs` to indices
        # into a new, all-unique array.
        # ----------------------------------------------------------
        l = sortedFaces.swapaxes(0,1).tolist() # arrange faces in pairs
        faceStrToIndex = {}
        facesIdxMap    = np.empty(len(l)).astype(int)
        numUniqueFaces = 0

        for i in range(len(l)):
            keyStr = ','.join([str(int(x)) for x in l[i]])
            if faceStrToIndex.has_key(keyStr): # duplicate
                facesIdxMap[i] = faceStrToIndex[keyStr] # refer to idx
            else: # not a duplicate
                faceStrToIndex[keyStr] = numUniqueFaces # establish idx
                facesIdxMap[i] = numUniqueFaces
                numUniqueFaces += 1
        # ----------------------------------------------------------

        # build all-unique numpy array
        noDups = np.empty((numUniqueFaces, self.dimension)).astype(int)
        for k in faceStrToIndex.keys():
            aFace = [int(x) for x in k.split(',')]
            noDups[faceStrToIndex[k]] = np.array(aFace)

        # finally, we build the cells from noDups-face indices
        numFacesPerCell = len(facesIdxMap[::numCells])
        cellsToFaces = np.empty((numCells, numFacesPerCell)).astype(int)
        for i in range(numCells):
            cellsToFaces[i] = facesIdxMap[i::numCells]

        # orient faces, cells for `Mesh` constructor
        return noDups.swapaxes(0,1), cellsToFaces.swapaxes(0,1)[::-1]

class GmshImporter2D(mesh2D.Mesh2D):
    """
    """
    def __init__(self, mshData, coordDimensions=2):
        mshFile = MshFile(mshData, 2)
        verts   = mshFile.vertexCoords
        faces   = mshFile.faces
        cells   = mshFile.cells
        print verts
        print faces
        print cells
        mesh2D.Mesh2D.__init__(self, vertexCoords=verts,
                                     faceVertexIDs=faces,
                                     cellFaceIDs=cells)
    def getCellVolumes(self):
        return abs(mesh.Mesh.getCellVolumes(self))
    
