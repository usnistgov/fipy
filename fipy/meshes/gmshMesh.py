#!/usr/bin/env python

## -*-Pyth-*-
# ###################################################################
#  FiPy - a finite volume PDE solver in Python
#
#  FILE: "gmshMesh.py"
#
#  Author: James O'Beirne <james.obeirne@nist.gov>
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
#  See the file "license.terms" for information on usage and
#  redistribution
#  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#
# ###################################################################
##

__docformat__ = 'restructuredtext'

import os
from subprocess import Popen, PIPE
import sys
import tempfile
from textwrap import dedent
import warnings
from distutils.version import StrictVersion

from fipy.tools import numerix as nx
from fipy.tools import parallelComm
from fipy.tools import serialComm
from fipy.tests.doctestPlus import register_skipper

from fipy.meshes.mesh import Mesh
from fipy.meshes.mesh2D import Mesh2D
from fipy.meshes.topologies.meshTopology import _MeshTopology

from fipy.tools.debug import PRINT

__all__ = ["openMSHFile", "openPOSFile",
           "Gmsh2D", "Gmsh2DIn3DSpace", "Gmsh3D",
           "GmshGrid2D", "GmshGrid3D"]

DEBUG = False

def _checkForGmsh():
    hasGmsh = True
    try:
        version = _gmshVersion(communicator=parallelComm)
        hasGmsh = version >= StrictVersion("2.0")
    except Exception:
        hasGmsh = False
    return hasGmsh

register_skipper(flag="GMSH",
                 test=_checkForGmsh,
                 why="`gmsh` cannot be found on the $PATH")

def parprint(str):
    if DEBUG:
        if parallelComm.procID == 0:
            print >> sys.stderr, str

class GmshException(Exception):
    pass

class MeshExportError(GmshException):
    pass

def gmshVersion(communicator=parallelComm):
    """Determine the version of Gmsh.

    We can't trust the generated msh file for the correct version number, so
    we have to retrieve it from the gmsh binary.
    """
    if communicator.procID == 0:
        while True:
            try:
                p = Popen(["gmsh", "--version"], stderr=PIPE)
            except OSError, e:
                verStr = None
                break

            try:
                out, verStr = p.communicate()
                verStr = verStr.decode('ascii').strip('\n')
                break
            except IOError:
                # some weird conflict with things like PyQT can cause
                # this to fail sometimes.
                # See http://thread.gmane.org/gmane.comp.python.enthought.devel/29362
                pass
    else:
        verStr = None

    return communicator.bcast(verStr)

def _gmshVersion(communicator=parallelComm):
    version = gmshVersion(communicator) or "0.0"
    return StrictVersion(version)

def openMSHFile(name, dimensions=None, coordDimensions=None, communicator=parallelComm, order=1, mode='r', background=None):
    """Open a Gmsh MSH file

    :Parameters:
      - `filename`: a string indicating gmsh output file
      - `dimensions`: an integer indicating dimension of mesh
      - `coordDimensions`: an integer indicating dimension of shapes
      - `order`: ???
      - `mode`: a string beginning with 'r' for reading and 'w' for writing.
        The file will be created if it doesn't exist when opened for writing;
        it will be truncated when opened for writing.
        Add a 'b' to the mode for binary files.
      - `background`: a `CellVariable` that specifies the desired characteristic
        lengths of the mesh cells
    """

    if order > 1:
        communicator = serialComm

    # Enforce gmsh version to be either >= 2 or 2.5, based on Nproc.
    version = _gmshVersion(communicator=communicator)
    if version < StrictVersion("2.0"):
        raise EnvironmentError("Gmsh version must be >= 2.0.")

    # If we're being passed a .msh file, leave it be. Otherwise,
    # we've gotta compile a .msh file from either (i) a .geo file,
    # or (ii) a gmsh script passed as a string.

    fileIsTemporary = False

    if mode.startswith('r'):
        if not os.path.exists(name):
            # we must have been passed a Gmsh script
            if communicator.procID == 0:
                (f, geoFile) = tempfile.mkstemp('.geo')
                file = os.fdopen(f, 'w')
                file.writelines(name)
                file.close()
            else:
                geoFile = None
            communicator.Barrier()
            geoFile = communicator.bcast(geoFile)
        else:
            # Gmsh isn't picky about file extensions,
            # so we peek at the start of the file to deduce the type
            f = open(name, 'r')
            filetype = f.readline().strip()
            f.close()
            if filetype == "$MeshFormat":
                geoFile = None
                mshFile = name
                gmshOutput = ""
            elif filetype == "$NOD":
                raise SyntaxError("Gmsh MSH file format version 1.0 is not supported")
            elif filetype == "$PostFormat":
                raise SyntaxError("Gmsh POS post-processing format cannot be used to generate a Mesh")
            else:
                # must be a Gmsh script file
                geoFile = name

        if geoFile is not None:
            gmshFlags = ["-%d" % dimensions, "-nopopup"]

            if communicator.Nproc > 1:
                if version < StrictVersion("2.5"):
                    warnstr = "Cannot partition with Gmsh version < 2.5. " \
                               + "Reverting to serial."
                    warnings.warn(warnstr, RuntimeWarning, stacklevel=2)
                    communicator = serialComm

                    dimensions = dimensions or coordDimensions

                    if dimensions is None:
                        raise ValueError("'dimensions' must be specified to generate a mesh from a geometry script")
                else: # gmsh version is adequate for partitioning
                    gmshFlags += ["-part", "%d" % communicator.Nproc]

            gmshFlags += ["-format", "msh"]

            if background is not None:
                if communicator.procID == 0:
                    f, bgmf = tempfile.mkstemp(suffix=".pos")
                    os.close(f)
                else:
                    bgmf = None
                bgmf = communicator.bcast(bgmf)
                f = openPOSFile(name=bgmf, mode='w')
                f.write(background)
                f.close()

                gmshFlags += ["-bgm", bgmf]

            if communicator.procID == 0:
                (f, mshFile) = tempfile.mkstemp('.msh')
                os.close(f)
                fileIsTemporary = True

                while True:
                    p = Popen(["gmsh", geoFile] + gmshFlags + ["-o", mshFile],
                              stdout=PIPE)

                    try:
                        gmshOutput, gmshError = p.communicate()
                        break
                    except IOError:
                        # some weird conflict with things like PyQT can cause
                        # this to fail sometimes.
                        # See http://thread.gmane.org/gmane.comp.python.enthought.devel/29362
                        pass

                if background is not None:
                    os.unlink(bgmf)

                if not os.path.exists(name):
                    os.unlink(geoFile)

                gmshOutput = gmshOutput.decode('ascii')

                parprint("gmsh out: %s" % gmshOutput)
            else:
                mshFile = None
                gmshOutput = ""

            mshFile = communicator.bcast(mshFile)
            gmshOutput = communicator.bcast(gmshOutput)
    elif mode.startswith('w'):
        mshFile = name
        gmshOutput = ""
    else:
        raise ValueError("mode string must begin with one of 'r' or 'w', not '%s'" % mode[0])

    return MSHFile(filename=mshFile,
                   dimensions=dimensions,
                   coordDimensions=coordDimensions,
                   communicator=communicator,
                   gmshOutput=gmshOutput,
                   mode=mode,
                   fileIsTemporary=fileIsTemporary)

def openPOSFile(name, communicator=parallelComm, mode='w'):
    """Open a Gmsh POS post-processing file
    """
    if not mode.startswith('w'):
        raise ValueError("mode string must begin with 'w', not '%s'" % mode[0])

    return POSFile(filename=name,
                   communicator=communicator,
                   mode=mode)

class GmshFile:
    def __init__(self, filename, communicator, mode, fileIsTemporary=False):
        self.filename = filename
        self.communicator = communicator
        self.mode = mode
        self.fileIsTemporary = fileIsTemporary

        self.formatWritten = False

        # open the .msh file
        if (hasattr(self.filename, "name")
            and hasattr(self.filename, "read")
            and hasattr(self.filename, "write")):
            self.fileobj = self.filename
            self.filename = self.fileobj.name
        else:
            self.fileobj = open(self.filename, mode=mode)

    def _getElementType(self, vertices, dimensions):
        if (vertices == 3 and dimensions == 2):
            return 2 ## triangle
        elif (vertices == 4 and dimensions == 2):
            return 3 ## quadrangle
        elif (vertices == 4 and dimensions == 3):
            return 4 ## tetrahedron
        elif (vertices == 8 and dimensions == 3):
            return 5 ## hexahedron
        elif (vertices == 6 and dimensions == 3):
            return 6 ## prism
        elif (vertices == 5 and dimensions == 3):
            return 7 ## pyramid
        else:
            raise MeshExportError, "Element type unsupported by Gmsh"

    def _orderVertices(self, vertexCoords, vertices):
        coordinates = nx.take(vertexCoords, vertices, axis=1)
        centroid = nx.add.reduce(coordinates, axis=1) / coordinates.shape[1]
        coordinates = coordinates - centroid[..., nx.newaxis]

        # to prevent div by zero
        coordinates = nx.where(coordinates == 0, 1.e-10, coordinates)

        # angles go from -pi / 2 to 3*pi / 2
        angles = nx.arctan(coordinates[1] / coordinates[0]) \
                + nx.where(coordinates[0] < 0, nx.pi, 0)
        sortorder = nx.argsort(angles)
        return nx.take(vertices, sortorder)

    def write(self, obj, time=0.0, timeindex=0):
        pass

    def close(self):
        self.fileobj.close()

    def __del__(self):
        if self.fileIsTemporary:
            os.unlink(self.filename)

class POSFile(GmshFile):
    def write(self, obj, time=0.0, timeindex=0):
        if not self.formatWritten:
            self._writeMeshFormat()
            self.formatWritten = True

        from fipy.variables.cellVariable import CellVariable

        if not isinstance(obj, CellVariable):
            raise TypeError("Unable to write %s to a POS file" % type(obj))

        coords = obj.mesh.vertexCoords
        dimensions, numNodes = coords.shape

        self._writeValues(var=obj, dimensions=dimensions, time=time, timeindex=timeindex)

    def _writeMeshFormat(self):
        versionNumber = 1.4
        sizeOfDouble = 8
        lines = ["$PostFormat\n",
                 "%g 0 %d\n" % (versionNumber, sizeOfDouble),
                 "$EndPostFormat\n"]
        self.fileobj.writelines(lines)

    def _writeValues(self, var, dimensions, time=0.0, timeindex=0):
        self.fileobj.write("$View\n")

        name = var.name
        if len(name) == 0:
            name = var.__class__.__name__
        # view-name nb-time-steps
        self.fileobj.write("%s %d\n" % (name.replace(" ", "_"), 1))

        # what a silly format

        # initialize all counts to zero
        nb_scalar_points = nb_vector_points = nb_tensor_points = 0
        nb_scalar_lines = nb_vector_lines = nb_tensor_lines = 0
        nb_scalar_triangles = nb_vector_triangles = nb_tensor_triangles = 0
        nb_scalar_quadrangles = nb_vector_quadrangles = nb_tensor_quadrangles = 0
        nb_scalar_tetrahedra = nb_vector_tetrahedra = nb_tensor_tetrahedra = 0
        nb_scalar_hexahedra = nb_vector_hexahedra = nb_tensor_hexahedra = 0
        nb_scalar_prisms = nb_vector_prisms = nb_tensor_prisms = 0
        nb_scalar_pyramids = nb_vector_pyramids = nb_tensor_pyramids = 0
        nb_scalar_lines2 = nb_vector_lines2 = nb_tensor_lines2 = 0
        nb_scalar_triangles2 = nb_vector_triangles2 = nb_tensor_triangles2 = 0
        nb_scalar_quadrangles2 = nb_vector_quadrangles2 = nb_tensor_quadrangles2 = 0
        nb_scalar_tetrahedra2 = nb_vector_tetrahedra2 = nb_tensor_tetrahedra2 = 0
        nb_scalar_hexahedra2 = nb_vector_hexahedra2 = nb_tensor_hexahedra2 = 0
        nb_scalar_prisms2 = nb_vector_prisms2 = nb_tensor_prisms2 = 0
        nb_scalar_pyramids2 = nb_vector_pyramids2 = nb_tensor_pyramids2 = 0
        nb_text2d = nb_text2d_chars = nb_text3d = nb_text3d_chars = 0

        nonOverlapping = nx.zeros((var.mesh.numberOfCells,), dtype=bool)
        nonOverlapping[..., var.mesh._localNonOverlappingCellIDs] = True

        cellTopology = var.mesh._cellTopology
        t = var.mesh.topology._elementTopology

        if var.rank == 0:
            nb_scalar_triangles = self.communicator.sum((cellTopology == t["triangle"]) & nonOverlapping)
            nb_scalar_quadrangles = self.communicator.sum((cellTopology == t["quadrangle"]) & nonOverlapping)
            nb_scalar_tetrahedra = self.communicator.sum((cellTopology == t["tetrahedron"]) & nonOverlapping)
            nb_scalar_hexahedra = self.communicator.sum((cellTopology == t["hexahedron"]) & nonOverlapping)
            nb_scalar_prisms = self.communicator.sum((cellTopology == t["prism"]) & nonOverlapping)
            nb_scalar_pyramids = self.communicator.sum((cellTopology == t["pyramid"]) & nonOverlapping)
        elif var.rank == 1:
            nb_vector_triangles = self.communicator.sum((cellTopology == t["triangle"]) & nonOverlapping)
            nb_vector_quadrangles = self.communicator.sum((cellTopology == t["quadrangle"]) & nonOverlapping)
            nb_vector_tetrahedra = self.communicator.sum((cellTopology == t["tetrahedron"]) & nonOverlapping)
            nb_vector_hexahedra = self.communicator.sum((cellTopology == t["hexahedron"]) & nonOverlapping)
            nb_vector_prisms = self.communicator.sum((cellTopology == t["prism"]) & nonOverlapping)
            nb_vector_pyramids = self.communicator.sum((cellTopology == t["pyramid"]) & nonOverlapping)
        elif var.rank == 2:
            nb_tensor_triangles = self.communicator.sum((cellTopology == t["triangle"]) & nonOverlapping)
            nb_tensor_quadrangles = self.communicator.sum((cellTopology == t["quadrangle"]) & nonOverlapping)
            nb_tensor_tetrahedra = self.communicator.sum((cellTopology == t["tetrahedron"]) & nonOverlapping)
            nb_tensor_hexahedra = self.communicator.sum((cellTopology == t["hexahedron"]) & nonOverlapping)
            nb_tensor_prisms = self.communicator.sum((cellTopology == t["prism"]) & nonOverlapping)
            nb_tensor_pyramids = self.communicator.sum((cellTopology == t["pyramid"]) & nonOverlapping)
        else:
            raise ValueError("rank must be 0, 1, or 2")


        self.fileobj.write(dedent("""\
            %(nb_scalar_points)d %(nb_vector_points)d %(nb_tensor_points)d
            %(nb_scalar_lines)d %(nb_vector_lines)d %(nb_tensor_lines)d
            %(nb_scalar_triangles)d %(nb_vector_triangles)d %(nb_tensor_triangles)d
            %(nb_scalar_quadrangles)d %(nb_vector_quadrangles)d %(nb_tensor_quadrangles)d
            %(nb_scalar_tetrahedra)d %(nb_vector_tetrahedra)d %(nb_tensor_tetrahedra)d
            %(nb_scalar_hexahedra)d %(nb_vector_hexahedra)d %(nb_tensor_hexahedra)d
            %(nb_scalar_prisms)d %(nb_vector_prisms)d %(nb_tensor_prisms)d
            %(nb_scalar_pyramids)d %(nb_vector_pyramids)d %(nb_tensor_pyramids)d
            %(nb_scalar_lines2)d %(nb_vector_lines2)d %(nb_tensor_lines2)d
            %(nb_scalar_triangles2)d %(nb_vector_triangles2)d %(nb_tensor_triangles2)d
            %(nb_scalar_quadrangles2)d %(nb_vector_quadrangles2)d %(nb_tensor_quadrangles2)d
            %(nb_scalar_tetrahedra2)d %(nb_vector_tetrahedra2)d %(nb_tensor_tetrahedra2)d
            %(nb_scalar_hexahedra2)d %(nb_vector_hexahedra2)d %(nb_tensor_hexahedra2)d
            %(nb_scalar_prisms2)d %(nb_vector_prisms2)d %(nb_tensor_prisms2)d
            %(nb_scalar_pyramids2)d %(nb_vector_pyramids2)d %(nb_tensor_pyramids2)d
            %(nb_text2d)d %(nb_text2d_chars)d %(nb_text3d)d %(nb_text3d_chars)d
            """ % locals()))

        # time-step-values
        self.fileobj.write("%g\n" % time)

        value = var.value

        vertexCoords = var.mesh.vertexCoords
        cellVertexIDs = var.mesh._orderedCellVertexIDs

        for shape in ("triangle", "quadrangle", "tetrahedron", "hexahedron", "prism", "pyramid"):
            for proc in range(self.communicator.Nproc):
                # write the elements in sequence to avoid pulling everything
                # to one processor.
                if proc == self.communicator.procID:
                    for i in ((cellTopology == t[shape]) & nonOverlapping).nonzero()[0]:
                        nodes = cellVertexIDs[..., i]
                        self._writeNodesAndValues(vertexCoords=vertexCoords,
                                                  nodes=nodes.compressed(),
                                                  value=value[..., i])
                self.communicator.Barrier()
                # need to synchronize all processes at the end of the file
                # this is ridiculously difficult to achieve
                self.fileobj.seek(self.communicator.bcast(self.fileobj.tell(), root=proc))

        self.fileobj.write("$EndView\n")

    # lists of double precision numbers giving the node coordinates and the
    # values associated with the nodes of the nb-scalar-points scalar points,
    # nb-vector-points vector points, ..., for each of the time-step-values.
    # For example, vector-triangle-value is defined as:
    #
    #           coord1-node1 coord1-node2 coord1-node3
    #           coord2-node1 coord2-node2 coord2-node3
    #           coord3-node1 coord3-node2 coord3-node3
    #           comp1-node1-time1 comp2-node1-time1 comp3-node1-time1
    #           comp1-node2-time1 comp2-node2-time1 comp3-node2-time1
    #           comp1-node3-time1 comp2-node3-time1 comp3-node3-time1
    #           comp1-node1-time2 comp2-node1-time2 comp3-node1-time2
    #           comp1-node2-time2 comp2-node2-time2 comp3-node2-time2
    #           comp1-node3-time2 comp2-node3-time2 comp3-node3-time2

    def _writeNodesAndValues(self, vertexCoords, nodes, value):
        # strip out masked values
        nodes = nodes[nodes != -1]
        numNodes = len(nodes)
        data = []
        dim = vertexCoords.shape[0]
        data = [[str(vertexCoords[..., j, node]) for node in nodes] for j in range(dim)]
        if dim == 2:
            data += [["0.0"] * numNodes]
        data += [[str(value)] * numNodes]
        self.fileobj.write("\n".join([" ".join(datum) for datum in data]) + "\n")

class MSHFile(GmshFile):
    """
    Class responsible for parsing a Gmsh file and then readying
    its contents for use by a `Mesh` constructor.

    Can handle a partitioned mesh based on `parallelComm.Nproc`. If partitioning,
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
                       communicator=parallelComm,
                       gmshOutput="",
                       mode='r',
                       fileIsTemporary=False):
        """
        :Parameters:
          - `filename`: a string indicating gmsh output file
          - `dimensions`: an integer indicating dimension of mesh
          - `coordDimensions`: an integer indicating dimension of shapes
          - `communicator`: ???
          - `gmshOutput`: output (if any) from Gmsh run that created .msh file
          - `mode`: a string beginning with 'r' for reading and 'w' for writing.
            The file will be created if it doesn't exist when opened for writing;
            it will be truncated when opened for writing.
            Add a 'b' to the mode for binary files.
          - `fileIsTemporary`: if `True`, `filename` should be cleaned up on deletion
        """
        self.dimensions = dimensions
        self.coordDimensions = coordDimensions
        self.gmshOutput = gmshOutput

        self.mesh = None
        self.meshWritten = False

        GmshFile.__init__(self, filename=filename, communicator=communicator, mode=mode, fileIsTemporary=fileIsTemporary)

    def _getMetaData(self):
        """
        Extracts gmshVersion, file-type, and data-size in that
        order.
        """
        self._seekForHeader("MeshFormat")
        metaData = self.fileobj.readline().split()
        self.fileobj.seek(0)
        return [float(x) for x in metaData]

    def _isolateData(self, title):
        """
        Gets all data between $[title] and $End[title], writes
        it out to its own file.
        """
        newF, newPath = tempfile.mkstemp(title, text=True)
        newF = os.fdopen(newF, 'w')
        try:
            self._seekForHeader(title)
        except Exception, e:
            newF.close()
            os.unlink(newPath)
            raise

        # extract the actual data within section
        while True:
            line = self.fileobj.readline()
            if ("$End%s" % title) not in line:
                newF.write(line)
            else: break

        self.fileobj.seek(0)
        newF.close()
        return newPath

    def _seekForHeader(self, title):
        """
        Iterate through a file until we end up at the section header
        for `title`. Function has obvious side-effects on `self.fileobj`.
        """
        while True:
            line = self.fileobj.readline()
            if len(line) == 0:
                raise EOFError("No `%s' header found!" % title)
                break
            elif (("$%s" % title) not in line):
                continue
            else:
                break # found header

    def _deriveCellsAndFaces(self, cellsToVertIDs, shapeTypes, numCells):
        """
        Uses element information obtained from `_parseElementFile` to deliver
        `facesToVertices` and `cellsToFaces`.
        """

        allShapes  = nx.unique(shapeTypes).tolist()
        maxFaces   = max([self.numFacesPerCell[x] for x in allShapes])

        # `cellsToFaces` must be padded with -1; see mesh.py
        currNumFaces = 0
        cellsToFaces = nx.ones((numCells, maxFaces), 'l') * -1
        facesDict    = {}
        uniqueFaces  = []

        # we now build `cellsToFaces` and `uniqueFaces`,
        # the latter will result in `facesToVertices`.
        for cellIdx in range(numCells):
            shapeType = shapeTypes[cellIdx]
            cell = cellsToVertIDs[cellIdx]

            if shapeType in [5, 12, 17]: # hexahedron
                faces = self._extractOrderedFaces(cell=cell,
                                                  faceOrderings=[[0, 1, 2, 3], # ordering of vertices gleaned from
                                                                 [4, 5, 6, 7], # a one-cube Grid3D example
                                                                 [0, 1, 5, 4],
                                                                 [3, 2, 6, 7],
                                                                 [0, 3, 7, 4],
                                                                 [1, 2, 6, 5]])
            elif shapeType in [6, 13, 18]: # prism
                faces = self._extractOrderedFaces(cell=cell,
                                                  faceOrderings=[[0, 1, 2],
                                                                 [5, 4, 3],
                                                                 [3, 4, 1, 0],
                                                                 [4, 5, 2, 1],
                                                                 [5, 3, 0, 2]])
            elif shapeType in [7, 14, 19]: # pyramid
                faces = self._extractOrderedFaces(cell=cell,
                                                  faceOrderings=[[0, 1, 2, 3],
                                                                 [0, 1, 4],
                                                                 [1, 2, 4],
                                                                 [2, 3, 4],
                                                                 [3, 0, 4]])
            else:
                if shapeType in [2, 9, 20, 21, 22, 23, 24, 25]:
                    faceLength = 2 # triangle
                elif shapeType in [3, 10, 16]:
                    faceLength = 2 # quadrangle
                elif shapeType in [4, 11, 29, 30, 31]:
                    faceLength = 3 # tetrahedron

                faces = self._extractRegularFaces(cell=cell,
                                                  faceLength=faceLength,
                                                  facesPerCell=self.numFacesPerCell[shapeType])

            for faceIdx in range(len(faces)):
                # NB: currFace is sorted for the key to spot duplicates
                currFace = faces[faceIdx]
                keyStr   = ' '.join([str(x) for x in sorted(currFace)])

                if keyStr in facesDict:
                    cellsToFaces[cellIdx][faceIdx] = facesDict[keyStr]
                else: # new face
                    facesDict[keyStr] = currNumFaces
                    cellsToFaces[cellIdx][faceIdx] = currNumFaces
                    uniqueFaces.append(currFace)
                    currNumFaces += 1

        # pad short faces with -1
        maxFaceLen = max([len(f) for f in uniqueFaces])
        uniqueFaces = [[-1] * (maxFaceLen - len(f)) + f for f in uniqueFaces]

        facesToVertices = nx.array(uniqueFaces, dtype=nx.INT_DTYPE)

        return facesToVertices.swapaxes(0,1)[::-1], cellsToFaces.swapaxes(0,1).copy('C'), facesDict

    def _translateNodesToVertices(self, entitiesNodes, vertexMap):
        """Translates entitiesNodes from Gmsh node IDs to `vertexCoords` indices.
        """
        entitiesVertices = []

        for entity in entitiesNodes:
            try:
                vertIndices = vertexMap[nx.array(entity)]
            except IndexError:
                vertIndices = nx.ones((len(entity),), 'l') * -1
            entitiesVertices.append(vertIndices)

        return entitiesVertices

    def _extractRegularFaces(self, cell, faceLength, facesPerCell):
        """Return faces for a regular poly(gon|hedron)

        Given `cell`, defined by node IDs, returns an array of
        `facesPerCell` faces of length `faceLength` in terms of nodes.
        """
        faces = []
        for i in range(facesPerCell):
            aFace = []
            for j in range(faceLength):
                aVertex = (i + j) % len(cell) # we may wrap
                aFace.append(int(cell[aVertex]))
            faces.append(aFace)
        return faces

    def _extractOrderedFaces(self, cell, faceOrderings):
        """Return faces for a 8-node hexahedron cell.
        """
        def orderingToFace(vertList):
            aFace = []
            for i in vertList:
                aFace.append(int(cell[i]))
            return aFace

        return [orderingToFace(o) for o in faceOrderings]

    def read(self):
        """
        0. Build cellsToVertices
        1. Recover needed vertexCoords and mapping from file using
           cellsToVertices
        2. Build cellsToVertIDs proper from vertexCoords and vertex map
        3. Build faces
        4. Build cellsToFaces

        Isolate relevant data into three files, store in
        `self.nodesPath` for $Nodes,
        `self.elemsPath` for $Elements.
        `self.namesFile` for $PhysicalNames.

        Returns vertexCoords, facesToVertexID, cellsToFaceID,
                cellGlobalIDMap, ghostCellGlobalIDMap.
        """
        self.version, self.fileType, self.dataSize = self._getMetaData()
        self.nodesPath = self._isolateData("Nodes")
        self.elemsPath = self._isolateData("Elements")
        try:
            self.namesPath = self._isolateData("PhysicalNames")
        except EOFError, e:
            self.namesPath = None

        try:
            if self.dimensions is None:
                nodesFile = open(self.nodesPath, 'r')
                nodesFile.readline() # skip number of nodes

                # We assume we have a 2D file unless we find a node
                # with a non-zero Z coordinate
                self.dimensions = 2
                for node in nodesFile:
                    line   = node.split()

                    newVert = [float(x) for x in line]

                    if newVert[2] != 0.0:
                        self.dimensions = 3
                        break

                nodesFile.close()

            self.coordDimensions = self.coordDimensions or self.dimensions

            # we need a conditional here so we don't pick up 2D shapes in 3D
            if self.dimensions == 2:
                self.numVertsPerFace = {1: 2, # 2-node line
                                        8: 2} # 3-node line
                self.numFacesPerCell = { 2: 3, # 3-node triangle (3 faces)
                                         9: 3, # 6-node triangle (we only read 1st 3)
                                        20: 3, # 9-node triangle (we only read 1st 3)
                                        21: 3, # 10-node triangle (we only read 1st 3)
                                        22: 3, # 12-node triangle (we only read 1st 3)
                                        23: 3, # 15-node triangle (we only read 1st 3)
                                        24: 3, # 15-node triangle (we only read 1st 3)
                                        25: 3, # 21-node triangle (we only read 1st 3)
                                         3: 4, # 4-node quadrangle (4 faces)
                                        10: 4, # 9-node quadrangle (we only read 1st 4)
                                        16: 4} # 8-node quadrangle (we only read 1st 4)
            elif self.dimensions == 3:
                self.numVertsPerFace = { 2: 3, # 3-node triangle (3 vertices)
                                         9: 3, # 6-node triangle (we only read 1st 3)
                                        20: 3, # 9-node triangle (we only read 1st 3)
                                        21: 3, # 10-node triangle (we only read 1st 3)
                                        22: 3, # 12-node triangle (we only read 1st 3)
                                        23: 3, # 15-node triangle (we only read 1st 3)
                                        24: 3, # 15-node triangle (we only read 1st 3)
                                        25: 3, # 21-node triangle (we only read 1st 3)
                                         3: 4, # 4-node quadrangle (4 vertices)
                                        10: 4, # 9-node quadrangle (we only read 1st 4)
                                        16: 4} # 8-node quadrangle (we only read 1st 4)
                self.numFacesPerCell = { 4: 4, # 4-node tetrahedron (4 faces)
                                        11: 4, # 10-node tetrahedron (we only read 1st 4)
                                        29: 4, # 20-node tetrahedron (we only read 1st 4)
                                        30: 4, # 35-node tetrahedron (we only read 1st 4)
                                        31: 4, # 56-node tetrahedron (we only read 1st 4)
                                         5: 6, # 8-node hexahedron (6 faces)
                                        12: 6, # 27-node tetrahedron (we only read 1st 6)
                                        17: 6, # 20-node tetrahedron (we only read 1st 6)
                                         6: 5, # 6-node prism (5 faces)
                                        13: 5, # 18-node prism (we only read 1st 6)
                                        18: 5, # 15-node prism (we only read 1st 6)
                                         7: 5, # 5-node pyramid (5 faces)
                                        14: 5, # 14-node pyramid (we only read 1st 5)
                                        19: 5} # 13-node pyramid (we only read 1st 5)
            else:
                raise GmshException("Mesh has fewer than 2 or more than 3 dimensions")

            parprint("Parsing elements.")
            (cellsData,
             ghostsData,
             facesData) = self._parseElementFile()

            cellsToGmshVerts = cellsData.nodes + ghostsData.nodes
            numCellsTotal    = len(cellsData.nodes) + len(ghostsData.nodes)
            allShapeTypes    = cellsData.shapes + ghostsData.shapes
            allShapeTypes    = nx.array(allShapeTypes)
            allShapeTypes    = nx.delete(allShapeTypes, nx.s_[numCellsTotal:])
            self.physicalCellMap = nx.array(cellsData.physicalEntities
                                            + ghostsData.physicalEntities)
            self.geometricalCellMap = nx.array(cellsData.geometricalEntities
                                               + ghostsData.geometricalEntities)

            if numCellsTotal < 1:
                errStr = "Gmsh hasn't produced any cells! Check your Gmsh code."
                errStr += "\n\nGmsh output:\n%s" % "".join(self.gmshOutput).rstrip()
                raise GmshException(errStr)

            parprint("Recovering coords.")
            parprint("numcells %d" % numCellsTotal)
            vertexCoords, vertIDtoIdx = self._vertexCoordsAndMap(cellsToGmshVerts)

            # translate Gmsh IDs to `vertexCoord` indices
            cellsToVertIDs = self._translateNodesToVertices(cellsToGmshVerts,
                                                            vertIDtoIdx)

            parprint("Building cells and faces.")
            (facesToV,
             cellsToF,
             facesDict) = self._deriveCellsAndFaces(cellsToVertIDs,
                                                    allShapeTypes,
                                                    numCellsTotal)

            # cell entities were easy to record on parsing
            # but we don't use Gmsh faces, so we need to correlate the nodes
            # that make up the Gmsh faces with the vertex IDs of the FiPy faces
            # so that we can check if any are named
            faceEntitiesDict = dict()

            # translate Gmsh IDs to `vertexCoord` indices
            facesToVertIDs = self._translateNodesToVertices(facesData.nodes,
                                                            vertIDtoIdx)

            for face, physicalEntity, geometricalEntity in zip(facesToVertIDs,
                                                               facesData.physicalEntities,
                                                               facesData.geometricalEntities):
                faceEntitiesDict[' '.join([str(x) for x in sorted(face)])] = (physicalEntity, geometricalEntity)

            self.physicalFaceMap = nx.zeros(facesToV.shape[-1:], 'l')
            self.geometricalFaceMap = nx.zeros(facesToV.shape[-1:], 'l')
            for face in facesDict.keys():
                # not all faces are necessarily tagged
                if face in faceEntitiesDict:
                    self.physicalFaceMap[facesDict[face]] = faceEntitiesDict[face][0]
                    self.geometricalFaceMap[facesDict[face]] = faceEntitiesDict[face][1]

            self.physicalNames = self._parseNamesFile()

        finally:
            os.unlink(self.nodesPath)
            os.unlink(self.elemsPath)
            if self.namesPath is not None:
                os.unlink(self.namesPath)

        # convert lists of cell vertices to a properly oriented masked array
        maxVerts = max([len(v) for v in cellsToVertIDs])
        # ticket:539 - NumPy 1.7 casts to array before concatenation and empty array defaults to float
        cellsToVertIDs = [nx.concatenate((v, nx.array([-1] * (maxVerts-len(v)), dtype=nx.INT_DTYPE))) for v in cellsToVertIDs]
        cellsToVertIDs = nx.MA.masked_equal(cellsToVertIDs, value=-1).swapaxes(0,1)

        parprint("Done with cells and faces.")
        return (vertexCoords, facesToV, cellsToF,
                cellsData.idmap, ghostsData.idmap,
                cellsToVertIDs)

    def write(self, obj, time=0.0, timeindex=0):
        if not self.formatWritten:
            self._writeMeshFormat()
            self.formatWritten = True

        from fipy.meshes.abstractMesh import AbstractMesh
        from fipy.variables.cellVariable import CellVariable

        if isinstance(obj, AbstractMesh):
            mesh = obj
        elif isinstance(obj, CellVariable):
            mesh = obj.mesh
        else:
            raise TypeError("Unable to write %s to a MSH file" % type(obj))

        if self.mesh is None:
            self.mesh = mesh
        elif mesh is not self.mesh:
            # Section 9.1 of the Gmsh manual:
            #   http://geuz.org/gmsh/doc/texinfo/#MSH-ASCII-file-format
            # says "Sections can be repeated in the same file",
            # but this seems to only apply to $NodeData, $ElementData,
            # and $ElementNodeData. Even if Gmsh understood it,
            # it's not worth the bookkeeping to be able to write more thean one.
            raise ValueError("Only one Mesh can be written to a MSH file")

        coords = mesh.vertexCoords
        dimensions, numNodes = coords.shape

        if not self.meshWritten:
            self._writeNodes(coords=coords, dimensions=dimensions, numNodes=numNodes)
            self._writeElements(mesh=self.mesh, coords=coords, dimensions=dimensions, numNodes=numNodes)

            self.meshWritten = True

        if isinstance(obj, CellVariable):
            self._writeValues(var=obj, dimensions=dimensions, time=time, timeindex=timeindex)

    def _writeMeshFormat(self):
        versionNumber = 2.2
        sizeOfDouble = 8
        lines = ["$MeshFormat\n",
                 "%f 0 %d\n" % (versionNumber, sizeOfDouble),
                 "$EndMeshFormat\n"]
        self.fileobj.writelines(lines)

    def _writeNodes(self, coords, dimensions, numNodes):
        self.fileobj.write("$Nodes\n")
        self.fileobj.write(str(numNodes) + '\n')

        for i in range(numNodes):
            self.fileobj.write("%s %s %s " % (str(i + 1),
                                              str(coords[0, i]),
                                              str(coords[1, i])))
            if dimensions == 2:
                self.fileobj.write("0 \n")
            elif dimensions == 3:
                self.fileobj.write(str(coords[2, i]))
                self.fileobj.write (" \n")
            else:
                raise MeshExportError, "Mesh has fewer than 2 or more than 3 dimensions"

        self.fileobj.write("$EndNodes\n")

    def _writeElements(self, mesh, coords, dimensions, numNodes):
        self.fileobj.write("$Elements\n")

        faceVertexIDs = mesh.faceVertexIDs
        cellFaceIDs = mesh.cellFaceIDs
        numCells = cellFaceIDs.shape[1]
        self.fileobj.write(str(numCells) + '\n')

        for i in range(numCells):
            ## build the vertex list
            vertexList = []

            for faceNum in cellFaceIDs[..., i]:
                """For more complicated meshes, some cells may have fewer
                faces than others. If this is the case, ignore the
                '--' entries."""
                if faceNum is nx.MA.masked:
                    continue
                for vertexNum in faceVertexIDs[..., faceNum]:
                    if vertexNum not in vertexList and vertexNum is not nx.MA.masked:
                        vertexList.append(vertexNum)

            if dimensions == 2:
                vertexList = self._orderVertices(coords, vertexList)

            numVertices = len(vertexList)
            elementType = self._getElementType(numVertices, dimensions)
            self.fileobj.write("%s %s 0 " % (str(i + 1), str(elementType)))

            self.fileobj.write(" ".join([str(a + 1) for a in vertexList]) + "\n")

        self.fileobj.write("$EndElements\n")

    def _writeValues(self, var, dimensions, time=0.0, timeindex=0):
        self.fileobj.write("$ElementData\n")

        # string-tags
        # "By default the first string-tag is interpreted as the name of the
        # post-processing view."
        self.fileobj.writelines([s + "\n" for s in ["1", "\"%s\"" % var.name]])

        # real-tags
        # "By default the first real-tag is interpreted as a time value
        # associated with the dataset."
        self.fileobj.writelines([s + "\n" for s in ["1", str(time)]])

        # integer-tags
        # "By default the first integer-tag is interpreted as a time step index
        # (starting at 0), the second as the number of field components of the
        # data in the view (1, 3 or 9), the third as the number of entities (nodes
        # or elements) in the view, and the fourth as the partition index for the
        # view data (0 for no partition)."
        self.fileobj.writelines([s + "\n" for s in ["4",
                                                    str(timeindex),
                                                    str(3**var.rank),
                                                    str(var.mesh.numberOfCells),
                                                    str(0)]])

        for i in range(var.mesh.numberOfCells):
            self.fileobj.write(" ".join([str(s) for s in [i+1] + list(var[..., i].value.flat)]) + "\n")

        self.fileobj.write("$EndElementData\n")

    def _vertexCoordsAndMap(self, cellsToGmshVerts):
        """
        Returns `vertexCoords` and mapping from Gmsh ID to `vertexCoords`
        indices (same as in MSHFile).

        Unlike parent, doesn't use genfromtxt
        because we want to avoid loading the entire msh file into memory.
        """
        allVerts     = [v for c in cellsToGmshVerts for v in c] # flatten
        allVerts     = nx.unique(nx.array(allVerts, dtype=nx.INT_DTYPE)) # remove dups
        allVerts     = nx.sort(allVerts)
        maxVertIdx   = allVerts[-1] + 1 # add one to offset zero
        vertGIDtoIdx = nx.ones(maxVertIdx, 'l') * -1 # gmsh ID -> vertexCoords idx
        vertexCoords = nx.empty((len(allVerts), self.coordDimensions))
        nodeCount    = 0

        # establish map. This works because allVerts is a sorted set.
        vertGIDtoIdx[allVerts] = nx.arange(len(allVerts))

        nodesFile = open(self.nodesPath, 'r')
        nodesFile.readline() # skip number of nodes

        # now we walk through node file with a sorted unique list of vertices
        # in hand. When we encounter 0th element in `allVerts`, save it
        # to `vertexCoords` then pop its ID off `allVerts`.
        for node in nodesFile:
            line   = node.split()
            nodeID = int(line[0])

            if nodeID == allVerts[nodeCount]:
                newVert = [float(x) for x in line[1:self.coordDimensions+1]]
                vertexCoords[nodeCount,:] = nx.array(newVert)
                nodeCount += 1

            if len(allVerts) == nodeCount:
                break

        nodesFile.close()

        # transpose for FiPy
        transCoords = vertexCoords.swapaxes(0,1)
        return transCoords, vertGIDtoIdx

    def _parseElementFile(self):
        """
        Return three objects, the first for non-ghost cells, the second for
        ghost cells, and the third for faces.

        All nastiness concerning ghost cell
        calculation is consolidated here: if we were ever to need to CALCULATE
        GHOST CELLS OURSELVES, the only code we'd have to change is in here.
        """

        def _parseTags(offset, currLineInts):
            if offset == -1:
                # if first valid shape
                offset = currLineInts[0]

            numTags = currLineInts[2]
            tags    = currLineInts[3:(3+numTags)]

            # the partition tags for don't seem to always be present
            # and don't always make much sense when they are

            if len(tags) >= 2:
                physicalEntity = tags.pop(0)
                geometricalEntity = tags.pop(0)
            else:
                physicalEntity = geometricalEntity = -1

            return offset, tags, physicalEntity, geometricalEntity


        cellsData = _ElementData()
        ghostsData = _ElementData()
        facesData = _ElementData()

        cellOffset = -1 # this will be subtracted from gmsh ID to obtain global ID
        faceOffset = -1 # this will be subtracted from gmsh ID to obtain global ID
        pid = self.communicator.procID + 1

        elemsFile = open(self.elemsPath, 'r')

        elemsFile.readline() # skip number of elements
        for el in elemsFile:
            currLineInts = [int(x) for x in el.split()]
            elemType     = currLineInts[1]

            if elemType in self.numFacesPerCell.keys():
                # element is a cell

                (cellOffset,
                 tags,
                 physicalEntity,
                 geometricalEntity) = _parseTags(offset=cellOffset,
                                                 currLineInts=currLineInts)
                currLineInts[0] -= cellOffset

                if len(tags) > 0:
                    # next item is a count
                    if tags[0] != len(tags) - 1:
                        warnings.warn("Partition count %d does not agree with number of remaining tags %d." % (tags[0], len(tags) - 1),
                                      SyntaxWarning, stacklevel=2)
                    tags = tags[1:]

                if self.communicator.Nproc > 1:
                    for tag in tags:
                        if -tag == pid:
                            # if we're collecting ghost cells and this is our ghost cell
                            ghostsData.add(currLine=currLineInts, elType=elemType,
                                           physicalEntity=physicalEntity,
                                           geometricalEntity=geometricalEntity)
                        elif tag == pid:
                            # el is in this processor's partition or we collect all cells
                            cellsData.add(currLine=currLineInts, elType=elemType,
                                          physicalEntity=physicalEntity,
                                          geometricalEntity=geometricalEntity)
                else:
                    # we collect all cells
                    cellsData.add(currLine=currLineInts, elType=elemType,
                                  physicalEntity=physicalEntity,
                                  geometricalEntity=geometricalEntity)
            elif elemType in self.numVertsPerFace.keys():
                # element is a face

                (faceOffset,
                 tags,
                 physicalEntity,
                 geometricalEntity) = _parseTags(offset=faceOffset,
                                                 currLineInts=currLineInts)
                currLineInts[0] -= faceOffset

                facesData.add(currLine=currLineInts, elType=elemType,
                              physicalEntity=physicalEntity,
                              geometricalEntity=geometricalEntity)

        elemsFile.close()

        return cellsData, ghostsData, facesData


    def _parseNamesFile(self):
        physicalNames = {
            0: dict(),
            1: dict(),
            2: dict(),
            3: dict()
        }
        if self.namesPath is not None:
            namesFile = open(self.namesPath, 'r')

            namesFile.readline() # skip number of elements
            for nm in namesFile:
                nm = nm.split()
                if self.version > 2.0:
                    dim = [int(nm.pop(0))]
                else:
                    # Gmsh format prior to 2.1 did not unambiguously tie
                    # physical names to physical entities of different dimensions
                    # http://article.gmane.org/gmane.comp.cad.gmsh.general/1601
                    dim = [0, 1, 2, 3]
                num = int(nm.pop(0))
                name = " ".join(nm)[1:-1]
                for d in dim:
                    physicalNames[d][name] = int(num)

            namesFile.close()

        return physicalNames

    def makeMapVariables(self, mesh):
        """Utility function to make MeshVariables that define different domains in the mesh
        """
        from fipy.variables.cellVariable import CellVariable
        from fipy.variables.faceVariable import FaceVariable

        self.physicalCellMap = CellVariable(mesh=mesh, value=self.physicalCellMap)
        self.geometricalCellMap = CellVariable(mesh=mesh, value=self.geometricalCellMap)
        self.physicalFaceMap = FaceVariable(mesh=mesh, value=self.physicalFaceMap)
        self.geometricalFaceMap = FaceVariable(mesh=mesh, value=self.geometricalFaceMap)

        physicalCells = dict()
        for name in self.physicalNames[self.dimensions].keys():
            physicalCells[name] = (self.physicalCellMap == self.physicalNames[self.dimensions][name])

        physicalFaces = dict()
        for name in self.physicalNames[self.dimensions-1].keys():
            physicalFaces[name] = (self.physicalFaceMap == self.physicalNames[self.dimensions-1][name])

        return (self.physicalCellMap,
                self.geometricalCellMap,
                physicalCells,
                self.physicalFaceMap,
                self.geometricalFaceMap,
                physicalFaces)

    def _test(self):
        """
        Test exporting

        >>> import os
        >>> import tempfile

        >>> from fipy import Grid2D, Tri2D, Grid3D, CylindricalGrid2D, CellVariable, doctest_raw_input
        >>> from fipy.meshes.uniformGrid2D import UniformGrid2D

        >>> dir = tempfile.mkdtemp()

        >>> g = Grid2D(nx = 10, ny = 10)
        >>> x, y = g.cellCenters
        >>> gvar = CellVariable(mesh=g, name="f(x,y)", value=1. / (x+y))
        >>> f = openMSHFile(name=os.path.join(dir, "g.msh"), mode='w') # doctest: +GMSH
        >>> f.write(g) # doctest: +GMSH
        >>> f.write(gvar) # doctest: +GMSH
        >>> f.close() # doctest: +GMSH

        >>> if __name__ == "__main__":
        ...     p = Popen(["gmsh", os.path.join(dir, "g.msh")]) # doctest: +GMSH
        ...     doctest_raw_input("Grid2D... Press enter.")

        >>> gg = GmshGrid2D(dx=1., dy=1., nx=10, ny=10) # doctest: +GMSH

        >>> f = openMSHFile(name=os.path.join(dir, "gg.msh"), mode='w') # doctest: +GMSH
        >>> f.write(gg) # doctest: +GMSH
        >>> f.close() # doctest: +GMSH

        >>> if __name__ == "__main__":
        ...     p = Popen(["gmsh", os.path.join(dir, "gg.msh")]) # doctest: +GMSH
        ...     doctest_raw_input("GmshGrid2D... Press enter.")

        >>> ug = UniformGrid2D(nx = 10, ny = 10)
        >>> f = openMSHFile(name=os.path.join(dir, "ug.msh"), mode='w') # doctest: +GMSH
        >>> f.write(ug) # doctest: +GMSH
        >>> f.close() # doctest: +GMSH

        >>> t = Tri2D(nx = 10, ny = 10)
        >>> f = openMSHFile(name=os.path.join(dir, "t.msh"), mode='w') # doctest: +GMSH
        >>> f.write(t) # doctest: +GMSH
        >>> f.close() # doctest: +GMSH

        >>> concat = ug + (t + ([10], [0]))
        >>> f = openMSHFile(name=os.path.join(dir, "concat.msh"), mode='w') # doctest: +GMSH
        >>> f.write(concat) # doctest: +GMSH
        >>> f.close() # doctest: +GMSH

        >>> if __name__ == "__main__":
        ...     p = Popen(["gmsh", os.path.join(dir, "concat.msh")]) # doctest: +GMSH
        ...     doctest_raw_input("Tri2D + Grid2D... Press enter.")

        >>> g3d = Grid3D(nx=10, ny=10, nz=30)
        >>> f = openMSHFile(name=os.path.join(dir, "g3d.msh"), mode='w') # doctest: +GMSH
        >>> f.write(g3d) # doctest: +GMSH
        >>> f.close() # doctest: +GMSH

        >>> if __name__ == "__main__":
        ...     p = Popen(["gmsh", os.path.join(dir, "g3d.msh")]) # doctest: +GMSH
        ...     doctest_raw_input("Grid3D... Press enter.")

        >>> cyl = CylindricalGrid2D(nx=10, ny=10)
        >>> f = openMSHFile(name=os.path.join(dir, "cyl.msh"), mode='w') # doctest: +GMSH
        >>> f.write(cyl) # doctest: +GMSH
        >>> f.close() # doctest: +GMSH

        >>> if __name__ == "__main__":
        ...     p = Popen(["gmsh", os.path.join(dir, "cyl.msh")]) # doctest: +GMSH
        ...     doctest_raw_input("CylindricalGrid2D... Press enter.")

        >>> import shutil
        >>> shutil.rmtree(dir)
        """
        pass

class _ElementData(object):
    """
    Bookkeeping for cells. Declared as own class for generality.

    "nodes": A Python list of the vertices that make up this element
    "shapes": A shapeTypes Python list
    "idmap": A Python list which maps vertexCoords idx -> global ID
    "physicalEntities": A Python list of the Gmsh physical entities each element is in
    "geometricalEntities": A Python list of the Gmsh geometrical entities each element is in
    """
    def __init__(self):
        self.nodes = []
        self.shapes = []
        self.idmap = [] # vertexCoords idx -> gmsh ID (global ID)
        self.physicalEntities = []
        self.geometricalEntities = []

    def add(self, currLine, elType, physicalEntity, geometricalEntity):
        numTags = currLine[2]
        self.nodes.append(currLine[(numTags+3):])
        self.shapes.append(elType)
        self.idmap.append(currLine[0])
        self.physicalEntities.append(physicalEntity)
        self.geometricalEntities.append(geometricalEntity)

class _GmshTopology(_MeshTopology):

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
        return nx.array(self.mesh.cellGlobalIDs)

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
        return nx.array(self.mesh.cellGlobalIDs + self.mesh.gCellGlobalIDs)

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
        return nx.arange(len(self.mesh.cellGlobalIDs))

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
        return nx.arange(len(self.mesh.cellGlobalIDs)
                         + len(self.mesh.gCellGlobalIDs))




class Gmsh2D(Mesh2D):
    """Construct a 2D Mesh using Gmsh

    >>> radius = 5.
    >>> side = 4.
    >>> squaredCircle = Gmsh2D('''
    ... // A mesh consisting of a square inside a circle inside a circle
    ...
    ... // define the basic dimensions of the mesh
    ...
    ... cellSize = 1;
    ... radius = %(radius)g;
    ... side = %(side)g;
    ...
    ... // define the compass points of the inner circle
    ...
    ... Point(1) = {0, 0, 0, cellSize};
    ... Point(2) = {-radius, 0, 0, cellSize};
    ... Point(3) = {0, radius, 0, cellSize};
    ... Point(4) = {radius, 0, 0, cellSize};
    ... Point(5) = {0, -radius, 0, cellSize};
    ...
    ... // define the compass points of the outer circle
    ...
    ... Point(6) = {-2*radius, 0, 0, cellSize};
    ... Point(7) = {0, 2*radius, 0, cellSize};
    ... Point(8) = {2*radius, 0, 0, cellSize};
    ... Point(9) = {0, -2*radius, 0, cellSize};
    ...
    ... // define the corners of the square
    ...
    ... Point(10) = {side/2, side/2, 0, cellSize/2};
    ... Point(11) = {-side/2, side/2, 0, cellSize/2};
    ... Point(12) = {-side/2, -side/2, 0, cellSize/2};
    ... Point(13) = {side/2, -side/2, 0, cellSize/2};
    ...
    ... // define the inner circle
    ...
    ... Circle(1) = {2, 1, 3};
    ... Circle(2) = {3, 1, 4};
    ... Circle(3) = {4, 1, 5};
    ... Circle(4) = {5, 1, 2};
    ...
    ... // define the outer circle
    ...
    ... Circle(5) = {6, 1, 7};
    ... Circle(6) = {7, 1, 8};
    ... Circle(7) = {8, 1, 9};
    ... Circle(8) = {9, 1, 6};
    ...
    ... // define the square
    ...
    ... Line(9) = {10, 13};
    ... Line(10) = {13, 12};
    ... Line(11) = {12, 11};
    ... Line(12) = {11, 10};
    ...
    ... // define the three boundaries
    ...
    ... Line Loop(1) = {1, 2, 3, 4};
    ... Line Loop(2) = {5, 6, 7, 8};
    ... Line Loop(3) = {9, 10, 11, 12};
    ...
    ... // define the three domains
    ...
    ... Plane Surface(1) = {2, 1};
    ... Plane Surface(2) = {1, 3};
    ... Plane Surface(3) = {3};
    ...
    ... // label the three domains
    ...
    ... // attention: if you use any "Physical" labels, you *must* label
    ... // all elements that correspond to FiPy Cells (Physical Surace in 2D
    ... // and Physical Volume in 3D) or Gmsh will not include them and FiPy
    ... // will not be able to include them in the Mesh.
    ...
    ... // note: if you do not use any labels, all Cells will be included.
    ...
    ... Physical Surface("Outer") = {1};
    ... Physical Surface("Middle") = {2};
    ... Physical Surface("Inner") = {3};
    ...
    ... // label the "north-west" part of the exterior boundary
    ...
    ... // note: you only need to label the Face elements
    ... // (Physical Line in 2D and Physical Surface in 3D) that correspond
    ... // to boundaries you are interested in. FiPy does not need them to
    ... // construct the Mesh.
    ...
    ... Physical Line("NW") = {5};
    ... ''' % locals()) # doctest: +GMSH

    It can be easier to specify certain domains and boundaries within Gmsh
    than it is to define the same domains and boundaries with FiPy expressions.

    Here we compare obtaining the same Cells and Faces using FiPy's
    parametric descriptions and Gmsh's labels.

    >>> x, y = squaredCircle.cellCenters # doctest: +GMSH

    >>> middle = ((x**2 + y**2 <= radius**2)
    ...           & ~((x > -side/2) & (x < side/2)
    ...               & (y > -side/2) & (y < side/2))) # doctest: +GMSH

    >>> print (middle == squaredCircle.physicalCells["Middle"]).all() # doctest: +GMSH
    True

    >>> X, Y = squaredCircle.faceCenters # doctest: +GMSH

    >>> NW = ((X**2 + Y**2 > (1.99*radius)**2)
    ...       & (X**2 + Y**2 < (2.01*radius)**2)
    ...       & (X <= 0) & (Y >= 0)) # doctest: +GMSH

    >>> print (NW == squaredCircle.physicalFaces["NW"]).all() # doctest: +GMSH
    True

    It is possible to direct Gmsh to give the mesh different densities in
    different locations

    >>> geo = '''
    ... // A mesh consisting of a square
    ...
    ... // define the corners of the square
    ...
    ... Point(1) = {1, 1, 0, 1};
    ... Point(2) = {0, 1, 0, 1};
    ... Point(3) = {0, 0, 0, 1};
    ... Point(4) = {1, 0, 0, 1};
    ...
    ... // define the square
    ...
    ... Line(1) = {1, 2};
    ... Line(2) = {2, 3};
    ... Line(3) = {3, 4};
    ... Line(4) = {4, 1};
    ...
    ... // define the boundary
    ...
    ... Line Loop(1) = {1, 2, 3, 4};
    ...
    ... // define the domain
    ...
    ... Plane Surface(1) = {1};
    ... '''

    >>> from fipy import CellVariable, numerix

    >>> std = []
    >>> bkg = None
    >>> for refine in range(4):
    ...     square = Gmsh2D(geo, background=bkg) # doctest: +GMSH
    ...     x, y = square.cellCenters # doctest: +GMSH
    ...     bkg = CellVariable(mesh=square, value=abs(x / 4) + 0.01) # doctest: +GMSH
    ...     std.append(numerix.std(numerix.sqrt(2 * square.cellVolumes) / bkg)) # doctest: +GMSH

    Check that the mesh is monotonically approaching the desired density

    >>> print numerix.greater(std[:-1], std[1:]).all() # doctest: +GMSH
    True

    and that the final density is close enough to the desired density

    >>> print std[-1] < 0.2 # doctest: +GMSH
    True

    The initial mesh doesn't have to be from Gmsh

    >>> from fipy import Tri2D

    >>> trisquare = Tri2D(nx=1, ny=1)
    >>> x, y = trisquare.cellCenters
    >>> bkg = CellVariable(mesh=trisquare, value=abs(x / 4) + 0.01)
    >>> std1 = numerix.std(numerix.sqrt(2 * trisquare.cellVolumes) / bkg)

    >>> square = Gmsh2D(geo, background=bkg) # doctest: +GMSH
    >>> x, y = square.cellCenters # doctest: +GMSH
    >>> bkg = CellVariable(mesh=square, value=abs(x / 4) + 0.01) # doctest: +GMSH
    >>> std2 = numerix.std(numerix.sqrt(2 * square.cellVolumes) / bkg) # doctest: +GMSH

    >>> print std1 > std2 # doctest: +GMSH
    True

    :Parameters:
      - `arg`: a string giving (i) the path to an MSH file, (ii) a path to a
         Gmsh geometry (".geo") file, or (iii) a Gmsh geometry script
      - `coordDimensions`: an integer indicating dimension of shapes
      - `order`: ???
      - `background`: a `CellVariable` that specifies the desired characteristic
        lengths of the mesh cells
    """

    def __init__(self,
                 arg,
                 coordDimensions=2,
                 communicator=parallelComm,
                 order=1,
                 background=None):

        self.mshFile = openMSHFile(arg,
                                   dimensions=2,
                                   coordDimensions=coordDimensions,
                                   communicator=communicator,
                                   order=order,
                                   mode='r',
                                   background=background)

        (verts,
         faces,
         cells,
         self.cellGlobalIDs,
         self.gCellGlobalIDs,
         self._orderedCellVertexIDs_data) = self.mshFile.read()

        self.mshFile.close()

        if communicator.Nproc > 1:
            self.globalNumberOfCells = communicator.sum(len(self.cellGlobalIDs))
            parprint("  I'm solving with %d cells total." % self.globalNumberOfCells)
            parprint("  Got global number of cells")

        Mesh2D.__init__(self, vertexCoords=verts,
                              faceVertexIDs=faces,
                              cellFaceIDs=cells,
                              communicator=communicator,
                              _TopologyClass=_GmshTopology)

        (self.physicalCellMap,
         self.geometricalCellMap,
         self.physicalCells,
         self.physicalFaceMap,
         self.geometricalFaceMap,
         self.physicalFaces) = self.mshFile.makeMapVariables(mesh=self)

        del self.mshFile

        parprint("Exiting Gmsh2D")

    def __setstate__(self, state):
        super(Gmsh2D, self).__setstate__(state)
        self.cellGlobalIDs = list(nx.arange(self.cellFaceIDs.shape[-1]))
        self.gCellGlobalIDs = []
        self.communicator = serialComm
        self.mshFile = None

    def _test(self):
        """
        First, we'll test Gmsh2D on a small circle with triangular
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
        ... ''') # doctest: +GMSH

        >>> print circ.cellVolumes[0] > 0 # doctest: +GMSH
        True

        Now we'll test Gmsh2D again, but on a rectangle.

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
        ... ''') # doctest: +GMSH

        >>> print rect.cellVolumes[0] > 0 # doctest: +GMSH
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
        ... ''') # doctest: +GMSH

        >>> print circle.cellVolumes[0] > 0 # doctest: +GMSH
        True

        >>> from fipy.tools import dump
        >>> f, tmpfile = dump.write(circle) # doctest: +GMSH
        >>> pickle_circle = dump.read(tmpfile, f) # doctest: +GMSH

        >>> print (pickle_circle.cellVolumes == circle.cellVolumes).all()
        ... # doctest: +GMSH, +SERIAL
        True

        >>> print (pickle_circle._globalOverlappingCellIDs == circle._globalOverlappingCellIDs).all()
        ... # doctest: +GMSH, +SERIAL
        True

        >>> (ftmp, mshFile) = tempfile.mkstemp('.msh')
        >>> os.close(ftmp)
        >>> f = openMSHFile(name=mshFile, mode='w') # doctest: +GMSH
        >>> f.write(circle) # doctest: +GMSH
        >>> f.close() # doctest: +GMSH

        >>> from fipy import doctest_raw_input
        >>> if __name__ == "__main__":
        ...     p = Popen(["gmsh", mshFile]) # doctest: +GMSH
        ...     doctest_raw_input("Circle... Press enter.")

        >>> os.remove(mshFile)

        >>> cmd = "Point(1) = {0, 0, 0, 0.05};"

        >>> Gmsh2D(cmd) #doctest: +IGNORE_EXCEPTION_DETAIL, +GMSH
        Traceback (most recent call last):
            ...
        GmshException: Gmsh hasn't produced any cells! Check your Gmsh code.







        Load a mesh consisting of a triangle and square

        >>> (fmsh, mshFile) = tempfile.mkstemp('.msh')
        >>> f = os.fdopen(fmsh, 'w')

        We need to do a little fancy footwork to account for multiple processes

        >>> partitions = [(i+1) * (-1 * (i != 0) + 1 * (i == 0)) for i in range(parallelComm.Nproc)]
        >>> numtags = 2 + 1 + len(partitions)
        >>> partitions = " ".join([str(i) for i in [parallelComm.Nproc] + partitions])

        >>> output = f.write('''$MeshFormat
        ... 2.2 0 8
        ... $EndMeshFormat
        ... $Nodes
        ... 5
        ... 1 0.0 0.0 0.0
        ... 2 1.0 0.0 0.0
        ... 3 2.0 0.0 0.0
        ... 4 0.0 1.0 0.0
        ... 5 1.0 1.0 0.0
        ... $EndNodes
        ... $Elements
        ... 2
        ... 1 3 %(numtags)s 99 2 %(partitions)s 1 2 5 4
        ... 2 2 %(numtags)s 98 2 %(partitions)s 2 3 5
        ... $EndElements
        ... ''' % locals())
        >>> f.close()

        >>> sqrTri = Gmsh2D(mshFile) # doctest: +GMSH

        >>> os.remove(mshFile)

        >>> print nx.allclose(sqrTri.cellVolumes, [1., 0.5]) # doctest: +GMSH
        True


        Write square and triangle volumes out as a POS file

        >>> from fipy import CellVariable
        >>> vol = CellVariable(mesh=sqrTri, value=sqrTri.cellVolumes) # doctest: +GMSH

        >>> if parallelComm.procID == 0:
        ...     (ftmp, posFile) = tempfile.mkstemp('.pos')
        ...     os.close(ftmp)
        ... else:
        ...     posFile = None
        >>> posFile = parallelComm.bcast(posFile)
        >>> f = openPOSFile(posFile, mode='w') # doctest: +GMSH
        >>> f.write(vol) # doctest: +GMSH
        >>> f.close() # doctest: +GMSH

        >>> f = open(posFile, mode='r')
        >>> print "".join(f.readlines()) # doctest: +GMSH
        $PostFormat
        1.4 0 8
        $EndPostFormat
        $View
        CellVariable 1
        0 0 0
        0 0 0
        1 0 0
        1 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0
        0 0 0 0
        0
        1.0 2.0 1.0
        0.0 0.0 1.0
        0.0 0.0 0.0
        0.5 0.5 0.5
        0.0 1.0 1.0 0.0
        0.0 0.0 1.0 1.0
        0.0 0.0 0.0 0.0
        1.0 1.0 1.0 1.0
        $EndView
        <BLANKLINE>

        >>> f.close()
        >>> parallelComm.Barrier()

        >>> if parallelComm.procID == 0:
        ...     os.remove(posFile)




        Load a mesh with no tags

        >>> (fmsh, mshFile) = tempfile.mkstemp('.msh')
        >>> f = os.fdopen(fmsh, 'w')

        >>> output = f.write('''$MeshFormat
        ... 2.2 0 8
        ... $EndMeshFormat
        ... $Nodes
        ... 3
        ... 1 0 0 0
        ... 2 0 1 0
        ... 3 1 1 0
        ... $EndNodes
        ... $Elements
        ... 1
        ... 1 2 0 1 3 2
        ... $EndElements
        ... ''' % locals())
        >>> f.close()

        >>> noTag = Gmsh2D(mshFile) # doctest: +GMSH, +SERIAL

        >>> os.remove(mshFile)

        """

class Gmsh2DIn3DSpace(Gmsh2D):
    """Create a topologically 2D Mesh in 3D coordinates using Gmsh

    :Parameters:
      - `arg`: a string giving (i) the path to an MSH file, (ii) a path to a
         Gmsh geometry (".geo") file, or (iii) a Gmsh geometry script
      - `coordDimensions`: an integer indicating dimension of shapes
      - `order`: ???
      - `background`: a `CellVariable` that specifies the desired characteristic
        lengths of the mesh cells
    """
    def __init__(self, arg, communicator=parallelComm, order=1, background=None):
        Gmsh2D.__init__(self,
                        arg,
                        coordDimensions=3,
                        communicator=communicator,
                        order=order,
                        background=background)

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
        ... ''').extrude(extrudeFunc=lambda r: 1.1 * r) # doctest: +GMSH

        >>> print sphere.cellVolumes[0] > 0 # doctest: +GMSH
        True

        >>> from fipy.tools import dump
        >>> f, tmpfile = dump.write(sphere) # doctest: +GMSH
        >>> pickle_sphere = dump.read(tmpfile, f) # doctest: +GMSH

        >>> print (pickle_sphere.cellVolumes == sphere.cellVolumes).all()
        ... # doctest: +GMSH, +SERIAL
        True

        >>> print (pickle_sphere._globalOverlappingCellIDs == sphere._globalOverlappingCellIDs).all()
        ... # doctest: +GMSH, +SERIAL
        True
        """
        pass

class Gmsh3D(Mesh):
    """Create a 3D Mesh using Gmsh

    :Parameters:
      - `arg`: a string giving (i) the path to an MSH file, (ii) a path to a
         Gmsh geometry (".geo") file, or (iii) a Gmsh geometry script
      - `order`: ???
      - `background`: a `CellVariable` that specifies the desired characteristic
        lengths of the mesh cells
    """
    def __init__(self, arg, communicator=parallelComm, order=1, background=None):
        self.mshFile  = openMSHFile(arg,
                                    dimensions=3,
                                    communicator=communicator,
                                    order=order,
                                    mode='r',
                                    background=background)

        (verts,
         faces,
         cells,
         self.cellGlobalIDs,
         self.gCellGlobalIDs,
         self._orderedCellVertexIDs_data) = self.mshFile.read()

        self.mshFile.close()

        Mesh.__init__(self, vertexCoords=verts,
                            faceVertexIDs=faces,
                            cellFaceIDs=cells,
                            communicator=communicator,
                            _TopologyClass=_GmshTopology)

        if self.communicator.Nproc > 1:
            self.globalNumberOfCells = self.communicator.sum(len(self.cellGlobalIDs))

        (self.physicalCellMap,
         self.geometricalCellMap,
         self.physicalCells,
         self.physicalFaceMap,
         self.geometricalFaceMap,
         self.physicalFaces) = self.mshFile.makeMapVariables(mesh=self)

        del self.mshFile

    def __setstate__(self, state):
        super(Gmsh3D, self).__setstate__(state)
        self.cellGlobalIDs = list(nx.arange(self.cellFaceIDs.shape[-1]))
        self.gCellGlobalIDs = []
        self.communicator = serialComm
        self.mshFile = None

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
        ... ''') # doctest: +GMSH

        >>> print prism.cellVolumes[0] > 0 # doctest: +GMSH
        True

        >>> from fipy.tools import dump
        >>> f, tmpfile = dump.write(prism) # doctest: +GMSH
        >>> pickle_prism = dump.read(tmpfile, f) # doctest: +GMSH

        >>> print (pickle_prism.cellVolumes == prism.cellVolumes).all()
        ... # doctest: +GMSH, +SERIAL
        True

        >>> print (pickle_prism._globalOverlappingCellIDs == prism._globalOverlappingCellIDs).all()
        ... # doctest: +GMSH, +SERIAL
        True

        Load a mesh consisting of a tetrahedron, prism, and pyramid

        >>> (fmsh, mshFile) = tempfile.mkstemp('.msh')
        >>> f = os.fdopen(fmsh, 'w')

        We need to do a little fancy footwork to account for multiple processes

        >>> partitions = [(i+1) * (-1 * (i != 0) + 1 * (i == 0)) for i in range(parallelComm.Nproc)]
        >>> numtags = 2 + 1 + len(partitions)
        >>> partitions = " ".join([str(i) for i in [parallelComm.Nproc] + partitions])

        >>> output = f.write('''$MeshFormat
        ... 2.2 0 8
        ... $EndMeshFormat
        ... $Nodes
        ... 8
        ... 1 0.0 0.0 0.0
        ... 2 1.0 0.0 0.0
        ... 3 1.0 1.0 0.0
        ... 4 1.0 0.0 1.0
        ... 5 3.0 0.0 0.0
        ... 6 3.0 1.0 0.0
        ... 7 3.0 0.0 1.0
        ... 8 2.0 0.5 -1.0
        ... $EndNodes
        ... $Elements
        ... 3
        ... 1 4 %(numtags)s 99 2 %(partitions)s 1 3 4 2
        ... 2 6 %(numtags)s 98 2 %(partitions)s 2 3 4 5 6 7
        ... 3 7 %(numtags)s 97 2 %(partitions)s 2 3 6 5 8
        ... $EndElements
        ... ''' % locals())
        >>> f.close()

        >>> tetPriPyr = Gmsh3D(mshFile) # doctest: +GMSH

        >>> os.remove(mshFile)

        >>> print nx.allclose(tetPriPyr.cellVolumes, [1./6, 1., 2./3]) # doctest: +GMSH
        True

        Write tetrahedron, prism, and pyramid volumes out as a POS file

        >>> from fipy import CellVariable
        >>> vol = CellVariable(mesh=tetPriPyr, value=tetPriPyr.cellVolumes, name="volume") # doctest: +GMSH

        >>> if parallelComm.procID == 0:
        ...     (ftmp, posFile) = tempfile.mkstemp('.pos')
        ...     os.close(ftmp)
        ... else:
        ...     posFile = None
        >>> posFile = parallelComm.bcast(posFile)
        >>> f = openPOSFile(posFile, mode='w') # doctest: +GMSH
        >>> f.write(vol) # doctest: +GMSH
        >>> f.close() # doctest: +GMSH

        >>> f = open(posFile, mode='r') # doctest: +GMSH
        >>> l = f.readlines() # doctest: +GMSH
        >>> f.close() # doctest: +GMSH
        >>> print "".join(l[:5]) # doctest: +GMSH
        $PostFormat
        1.4 0 8
        $EndPostFormat
        $View
        volume 1
        <BLANKLINE>

        >>> print l[-1] # doctest: +GMSH
        $EndView
        <BLANKLINE>

        Py3k writes the numbers at a different precision

        >>> from fipy import numerix

        >>> a1 = numerix.fromstring("".join(l[5:-1]), sep=" ") # doctest: +GMSH
        >>> a2 = numerix.fromstring('''
        ...  0 0 0
        ...  0 0 0
        ...  0 0 0
        ...  0 0 0
        ...  1 0 0
        ...  0 0 0
        ...  1 0 0
        ...  1 0 0
        ...  0 0 0
        ...  0 0 0
        ...  0 0 0
        ...  0 0 0
        ...  0 0 0
        ...  0 0 0
        ...  0 0 0
        ...  0 0 0 0
        ...  0
        ...  0.0 1.0 1.0 1.0
        ...  0.0 1.0 0.0 0.0
        ...  0.0 0.0 1.0 0.0
        ...  0.16666666666666663 0.16666666666666663 0.16666666666666663 0.16666666666666663
        ...  1.0 1.0 1.0 3.0 3.0 3.0
        ...  0.0 1.0 0.0 0.0 1.0 0.0
        ...  0.0 0.0 1.0 0.0 0.0 1.0
        ...  1.0 1.0 1.0 1.0 1.0 1.0
        ...  1.0 1.0 3.0 3.0 2.0
        ...  0.0 1.0 1.0 0.0 0.5
        ...  0.0 0.0 0.0 0.0 -1.0
        ...  0.6666666666666666 0.6666666666666666 0.6666666666666666 0.6666666666666666 0.6666666666666666
        ...  ''', sep=" ")
        >>> print numerix.allclose(a1, a2) # doctest: +GMSH
        True

        >>> if parallelComm.procID == 0:
        ...     os.remove(posFile)
        """

class GmshGrid2D(Gmsh2D):
    """Should serve as a drop-in replacement for Grid2D."""
    def __init__(self, dx=1., dy=1., nx=1, ny=None,
                 coordDimensions=2, communicator=parallelComm, order=1):
        self.dx = dx
        self.dy = dy or dx
        self.nx = nx
        self.ny = ny or nx

        arg = self._makeGridGeo(self.dx, self.dy, self.nx, self.ny)

        Gmsh2D.__init__(self, arg, coordDimensions, communicator, order, background=None)

    @property
    def _meshSpacing(self):
        return nx.array((self.dx,self.dy))[...,nx.newaxis]

    def _makeGridGeo(self, dx, dy, nx, ny):
        height = ny * dy
        width  = nx * dx
        numLayers = int(ny / float(dy))

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

        >>> yogmsh = GmshGrid2D(dx=5, dy=5, nx=5, ny=5, communicator=serialComm) # doctest: +GMSH

        >>> yogrid = Grid2D(dx=5, dy=5, nx=5, ny=5, communicator=serialComm)

        >>> numerix.allclose(yogmsh._faceAreas, yogrid._faceAreas) # doctest: +GMSH
        True

        >>> yogmsh.cellCenters.value.size == yogrid.cellCenters.value.size # doctest: +GMSH
        True

        >>> mesh = GmshGrid2D(nx=2, ny=2) # doctest: +GMSH

        >>> mesh.numberOfCells == 4 # doctest: +GMSH
        True

        >>> len(mesh.faceCenters[0]) == 12 # doctest: +GMSH
        True
        """

class GmshGrid3D(Gmsh3D):
    """Should serve as a drop-in replacement for Grid3D."""
    def __init__(self, dx=1., dy=1., dz=1., nx=1, ny=None, nz=None,
                 communicator=parallelComm, order=1):
        self.dx = dx
        self.dy = dy or dx
        self.dz = dz or dx

        self.nx = nx
        self.ny = ny or nx
        self.nz = nz or nx

        arg = self._makeGridGeo(self.dx, self.dy, self.dz,
                                self.nx, self.ny, self.nz)

        Gmsh3D.__init__(self, arg, communicator=communicator, order=order)

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
        ...                     communicator=serialComm) # doctest: +GMSH

        >>> yogrid = Grid3D(dx=5, dy=5, dz=5, nx=5, ny=5, nz=5,
        ...                 communicator=serialComm)

        >>> yogmsh.cellCenters.value.size == yogrid.cellCenters.value.size # doctest: +GMSH
        True

        >>> numerix.allclose(yogmsh._faceAreas, yogrid._faceAreas) # doctest: +GMSH
        True

        >>> numerix.allclose(yogmsh._faceAreas, yogrid._faceAreas) # doctest: +GMSH
        True

        >>> mesh = GmshGrid3D(nx=2, ny=2, nz=2) # doctest: +GMSH

        >>> ccs = [[ 0.5,  0.5,  0.5,  0.5,  1.5,  1.5,  1.5,  1.5],
        ...    [ 0.5,  0.5,  1.5,  1.5,  0.5,  0.5,  1.5,  1.5],
        ...    [ 0.5,  1.5,  0.5,  1.5,  0.5,  1.5,  0.5,  1.5]]

        >>> len(mesh.cellCenters.value[0]) == 8 # doctest: +GMSH
        True

        >>> faceAreas = [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
        ...           1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
        ...           1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.]

        >>> nx.allclose(mesh._faceAreas, faceAreas) # doctest: +GMSH
        True

        >>> cellAreas = [[ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
        ...            [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
        ...            [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
        ...            [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
        ...            [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
        ...            [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.]]

        >>> nx.allclose(mesh._cellAreas, cellAreas) # doctest: +GMSH
        True
        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
