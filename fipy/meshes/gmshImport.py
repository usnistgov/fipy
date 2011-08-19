#!/usr/bin/env python
 
## -*-Pyth-*-
# ###################################################################
#  FiPy - a finite volume PDE solver in Python
#
#  FILE: "gmshImport.py"
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
import subprocess as subp
import sys
import tempfile
from textwrap import dedent
import warnings

from fipy.tools import numerix as nx
from fipy.tools import parallel
from fipy.tools import serial
from fipy.tools.decorators import getsetDeprecated
from fipy.meshes.mesh import Mesh
from fipy.meshes.mesh2D import Mesh2D

DEBUG = False


def parprint(str):
    if DEBUG:
        if parallel.procID == 0:
            print >> sys.stderr, str

class GmshException(Exception):
    pass
    
def _gmshVersion():
    """Determine the version of Gmsh.
    
    We can't trust the generated msh file for the correct version number, so
    we have to retrieve it from the gmsh binary.
    """
    import re

    verStr = "\n".join(subp.Popen("gmsh --version", 
                                  stderr=subp.PIPE, shell=True).stderr.readlines())

    m = re.search(r'\d+.\d+', verStr)

    if m:
        return float(m.group(0))
    else:
        return 0
    
def openMSHFile(name, dimensions=None, coordDimensions=None, communicator=parallel, order=1, mode='r', background=None):
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
    """
    
    if order > 1:
        communicator = serial

    # Enforce gmsh version to be either >= 2 or 2.5, based on Nproc.
    gmshVersion = _gmshVersion()
    if gmshVersion < 2.0:
        raise EnvironmentError("Gmsh version must be >= 2.0.")

    # If we're being passed a .msh file, leave it be. Otherwise,
    # we've gotta compile a .msh file from either (i) a .geo file, 
    # or (ii) a gmsh script passed as a string.
    
    if mode.startswith('r'):
        if not os.path.exists(name):
            # we must have been passed a Gmsh script
            (f, geoFile) = tempfile.mkstemp('.geo')
            file = open(geoFile, 'w')
            file.writelines(name)
            file.close() 
            os.close(f)
        else:
            # Gmsh isn't picky about file extensions, 
            # so we peek at the start of the file to deduce the type
            f = open(name, 'r')
            filetype = f.readline()
            f.close()
            if filetype == "$MeshFormat":
                geoFile = None
                mshFile = name
                gmshOutput = None
            elif filetype == "$NOD":
                raise SyntaxError("Gmsh MSH file format version 1.0 is not supported")
            elif filetype == "$PostFormat":
                raise SyntaxError("Gmsh POS post-processing format cannot be used to generate a Mesh")
            else:
                # must be a Gmsh script file
                geoFile = name
                
        if geoFile is not None:
            if communicator.Nproc > 1:
                if gmshVersion < 2.5:
                    warnstr = "Cannot partition with Gmsh version < 2.5. " \
                               + "Reverting to serial."
                    warnings.warn(warnstr, RuntimeWarning, stacklevel=2)
                    communicator = serial
                    
                    dimensions = dimensions or coordDimensions
                    
                    if dimensions is None:
                        raise ValueError("'dimensions' must be specified to generate a mesh from a geometry script")

                    gmshFlags = "-%d" % (dimensions)
                else: # gmsh version is adequate for partitioning
                    gmshFlags = "-%d -part %d" % (dimensions,
                                                  communicator.Nproc)
            else: # we're running serial
                gmshFlags = "-%d" % (dimensions)
            
            gmshFlags += " -format msh"
            
            if background is not None:
                bgmf = tempfile.NamedTemporaryFile(mode='w', suffix='.pos', delete=False)
                f = openPOSFile(name=bgmf, mode='w')
                f.write(background)
                f.close()

                gmshFlags += " -bgm %s" % bgmf.name

            (f, mshFile) = tempfile.mkstemp('.msh')
            gmshOutput = subp.Popen("gmsh %s %s -o %s"
                                    % (geoFile, gmshFlags, mshFile),
                                    stdout=subp.PIPE, shell=True).stdout.readlines()
            parprint("gmsh out: %s" % gmshOutput)
            os.close(f)
    elif mode.startswith('w'):
        mshFile = name
        gmshOutput = None
    else:
        raise ValueError("mode string must begin with one of 'r' or 'w', not '%s'" % mode[0])
                
    return MSHFile(filename=mshFile, 
                   dimensions=dimensions, 
                   coordDimensions=coordDimensions,
                   communicator=communicator, 
                   gmshOutput=gmshOutput,
                   mode=mode)
    
def openPOSFile(name, communicator=parallel, mode='w'):
    """Open a Gmsh POS post-processing file
    """
    if not mode.startswith('w'):
        raise ValueError("mode string must begin with 'w', not '%s'" % mode[0])
                
    return POSFile(filename=name, 
                   communicator=communicator, 
                   mode=mode)
    
class GmshFile:
    def __init__(self, filename, communicator, mode):
        self.filename = filename
        self.communicator = communicator
        self.mode = mode
        
        self.formatWritten = False
        
        # open the .msh file
        if (isinstance(self.filename, file) 
            or (hasattr(self.filename, "file") 
                and isinstance(self.filename.file, file))):
            self.fileobj = self.filename
            self.filename = self.fileobj.name
        else:
            self.fileobj = open(name=self.filename, mode=mode)

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
              
        # view-name nb-time-steps
        self.fileobj.write("%s %d\n" % (var.name.replace(" ", "_"), 1))
        
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
        
        maxFacesPerCell = var.mesh.cellFaceIDs.shape[0]
        if nx.MA.is_masked(var.mesh.cellFaceIDs):
            facesPerCell = (~nx.MA.getmask(var.mesh.cellFaceIDs)).sum(axis=0)
        else:
            facesPerCell = nx.empty((var.mesh.numberOfCells,), dtype=int)
            facesPerCell[:] = maxFacesPerCell
            
        cellFaceVertices = nx.take(var.mesh.faceVertexIDs, var.mesh.cellFaceIDs, axis=1)
        # argsort does not work as documented:
        #   Array of indices that sort a along the specified axis. In other words, a[index_array] yields a sorted a.
        # No, it's not.
        # We need both the sorted values and the sort order.
        if nx.MA.is_masked(cellFaceVertices):
            nodesPerFace = (~cellFaceVertices.mask).sum(axis=0)
        else:
            nodesPerFace = nx.empty(cellFaceVertices.shape[1:], dtype=int)
            nodesPerFace[:] = var.mesh.faceVertexIDs.shape[0]
        faceOrder = nx.argsort(nodesPerFace, axis=0)[::-1]
        nodesPerFace = nx.sort(nodesPerFace, axis=0)[::-1]
        
        if dimensions == 2:
            triangles = (facesPerCell == 3)
            quadrangles = (facesPerCell == 4)
            tetrahedra = nx.array([False])
            hexahedra = nx.array([False])
            prisms = nx.array([False])
            pyramids = nx.array([False])
            if var.rank == 0:
                nb_scalar_triangles = triangles.sum()
                nb_scalar_quadrangles = quadrangles.sum()
            elif var.rank == 1:
                nb_vector_triangles = triangles.sum()
                nb_vector_quadrangles = quadrangles.sum()
            elif var.rank == 2:
                nb_tensor_triangles = triangles.sum()
                nb_tensor_quadrangles = quadrangles.sum()
            else:
                raise ValueError("rank must be 0, 1, or 2")
        elif dimensions == 3:
            def faceNodeCount(counts):
                counts = counts + [0] * (maxFacesPerCell - len(counts))
                return nx.array(counts)[..., nx.newaxis]
                
            triangles = nx.array([False])
            quadrangles = nx.array([False])
            tetrahedra = ((facesPerCell == 4) 
                          & (nodesPerFace == faceNodeCount([3, 3, 3, 3])).any(axis=0))
            hexahedra = ((facesPerCell == 6) 
                         & (nodesPerFace == faceNodeCount([4, 4, 4, 4, 4, 4])).any(axis=0))
            prisms = ((facesPerCell == 5) 
                      & (nodesPerFace == faceNodeCount([4, 4, 4, 3, 3])).any(axis=0))
            pyramids = ((facesPerCell == 5) 
                        & (nodesPerFace == faceNodeCount([4, 3, 3, 3, 3])).any(axis=0))
        else:
            raise ValueError("dimensions must be 2 or 3")
            
        if var.rank == 0:
            nb_scalar_triangles = triangles.sum()
            nb_scalar_quadrangles = quadrangles.sum()
            nb_scalar_tetrahedra = tetrahedra.sum()
            nb_scalar_hexahedra = hexahedra.sum()
            nb_scalar_prisms = prisms.sum()
            nb_scalar_pyramids = pyramids.sum()
        elif var.rank == 1:
            nb_vector_triangles = triangles.sum()
            nb_vector_quadrangles = quadrangles.sum()
            nb_vector_tetrahedra = tetrahedra.sum()
            nb_vector_hexahedra = hexahedra.sum()
            nb_vector_prisms = prisms.sum()
            nb_vector_pyramids = pyramids.sum()
        elif var.rank == 2:
            nb_tensor_triangles = triangles.sum()
            nb_tensor_quadrangles = quadrangles.sum()
            nb_tensor_tetrahedra = tetrahedra.sum()
            nb_tensor_hexahedra = hexahedra.sum()
            nb_tensor_prisms = prisms.sum()
            nb_tensor_pyramids = pyramids.sum()
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

        for i in triangles.nonzero()[0]:
            triangle = cellFaceVertices[..., faceOrder[..., i], i]
            self._writeTriangle(mesh=var.mesh, triangle=triangle, value=value[..., i])
                                
        for i in quadrangles.nonzero()[0]:
            quadrangle = cellFaceVertices[..., faceOrder[..., i], i]
            self._writeQuadrangle(mesh=var.mesh, quadrangle=quadrangle, value=value[..., i])

        for i in tetrahedra.nonzero()[0]:
            tetrahedron = cellFaceVertices[..., faceOrder[..., i], i]
            self._writeTetrahedron(mesh=var.mesh, tetrahedron=tetrahedron, value=value[..., i])

        for i in hexahedra.nonzero()[0]:
            hexahedron = cellFaceVertices[..., faceOrder[..., i], i]
            self._writeHexahedron(mesh=var.mesh, hexahedron=hexahedron, value=value[..., i])

        for  i in prisms.nonzero()[0]:
            prism = cellFaceVertices[..., faceOrder[..., i], i]
            self._writePrism(mesh=var.mesh, prism=prism, value=value[..., i])

        for i in pyramids.nonzero()[0]:
            pyramid = cellFaceVertices[..., faceOrder[..., i], i]
            self._writePyramid(mesh=var.mesh, pyramid=pyramid, value=value[..., i])

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
    
    def _writeTriangle(self, mesh, triangle, value):
        # triangle is defined by one face and the remaining point 
        # from either of the other faces
        nodes = triangle[..., 0]
        nodes = nx.concatenate((nodes, triangle[~nx.in1d(triangle[..., 1], nodes), 1]))
                        
        self._writeNodesAndValues(mesh=mesh, nodes=nodes, value=value)

    def _writeQuadrangle(self, mesh, quadrangle, value):
        # quadrangle is defined by one face and the opposite face
        face0 = quadrangle[..., 0]
        for face1 in quadrangle[..., 1:].swapaxes(0, 1):
            if not nx.in1d(face1, face0).any():
                break
                
        # need to ensure face1 is oriented same way as face0
        cross01 = nx.cross((mesh.vertexCoords[..., face1[0]] 
                            - mesh.vertexCoords[..., face0[0]]),
                           (mesh.vertexCoords[..., face1[1]] 
                            - mesh.vertexCoords[..., face0[0]]))
        
        cross10 = nx.cross((mesh.vertexCoords[..., face0[0]] 
                            - mesh.vertexCoords[..., face1[0]]),
                           (mesh.vertexCoords[..., face0[1]] 
                            - mesh.vertexCoords[..., face1[0]]))

        if ((mesh.dim == 2 and cross01 * cross10 < 0)
            or (mesh.dim == 3 and nx.dot(cross01, cross10) < 0)):
            face1 = face1[::-1]
                
        nodes = nx.concatenate((face0, face1))
         
        self._writeNodesAndValues(mesh=mesh, nodes=nodes, value=value)
        
    def _writeTetrahedron(self, mesh, tetrahedron, value):
        # tetrahedron is defined by one face and the remaining point 
        # from any of the other faces
        nodes = tetrahedron[..., 0]
        nodes = nx.concatenate((nodes, tetrahedron[~nx.in1d(tetrahedron[..., 1], nodes), 1]))
                
        self._writeNodesAndValues(mesh=mesh, nodes=nodes, value=value)

    @profile
    def _reorientFace(self, mesh, face0, face1):
        # need to ensure face1 is oriented same way as face0
        cross01 = nx.cross((mesh.vertexCoords[..., face0[1]] 
                            - mesh.vertexCoords[..., face0[0]]),
                           (mesh.vertexCoords[..., face0[2]] 
                            - mesh.vertexCoords[..., face0[0]]))
                    
        cross10 = nx.cross((mesh.vertexCoords[..., face1[1]] 
                            - mesh.vertexCoords[..., face1[0]]),
                           (mesh.vertexCoords[..., face1[2]] 
                            - mesh.vertexCoords[..., face1[0]]))
                            
        if nx.dot(cross01, cross10) < 0:
            face1 = face1[::-1]
            
        return face1

    @profile
    def _writeHexahedron(self, mesh, hexahedron, value):
        # hexahedron is defined by one face and the opposite face
        face0 = hexahedron[..., 0]
        for face1 in hexahedron[..., 1:].swapaxes(0, 1):
            if not nx.in1d(face1, face0).any():
                break

        face1 = self._reorientFace(mesh=mesh, face0=face0, face1=face1)

        nodes = nx.concatenate((face0, face1))

        self._writeNodesAndValues(mesh=mesh, nodes=nodes, value=value)

    def _writePrism(self, mesh, prism, value):
        # prism is defined by the two three-sided faces 
        face0 = pyramid[..., 3]
        face1 = pyramid[..., 4]
        
        face1 = self._reorientFace(mesh=mesh, face0=face0, face1=face1)

        nodes = nx.concatenate((face0, face1))

        self._writeNodesAndValues(mesh=mesh, nodes=nodes, value=value)


    def _writePyramid(self, mesh, pyramid, value):
        # pyramid is defined by four-sided face and the remaining point 
        # from any of the other faces
        nodes = pyramid[..., 0]
        nodes = nx.concatenate((nodes, pyramid[~nx.in1d(pyramid[..., 1], nodes), 1]))
                       
        self._writeNodesAndValues(mesh=mesh, nodes=nodes, value=value)
        
    def _writeNodesAndValues(self, mesh, nodes, value):
        numNodes = len(nodes)
        data = []
        # this might be calculated (e.g. UniformGrid...)
        vertexCoords = mesh.vertexCoords
        data = [[str(vertexCoords[..., j, nodes[i]]) for i in range(numNodes)] for j in range(mesh.dim)]
        if mesh.dim == 2:
            data += [["0.0"] * numNodes]
        data += [[str(value)] * numNodes]
        self.fileobj.write("\n".join([" ".join(datum) for datum in data]) + "\n")



class MSHFile(GmshFile):
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
                       gmshOutput=None,
                       mode='r'):
        """
        :Parameters:
          - `filename`: a string indicating gmsh output file
          - `dimensions`: an integer indicating dimension of mesh
          - `coordDimensions`: an integer indicating dimension of shapes
          - `communictator`: ???
          - `gmshOutput`: output (if any) from Gmsh run that created .msh file
          - `mode`: a string beginning with 'r' for reading and 'w' for writing. 
            The file will be created if it doesn't exist when opened for writing; 
            it will be truncated when opened for writing.  
            Add a 'b' to the mode for binary files.
        """
        self.dimensions = dimensions
        self.coordDimensions = coordDimensions
        self.gmshOutput = gmshOutput
        
        self.mesh = None
        self.meshWritten = False

        GmshFile.__init__(self, filename=filename, communictator=communictator, mode=mode)
        
    def _getMetaData(self):
        """
        Extracts gmshVersion, file-type, and data-size in that
        order.
        """
        self._seekForHeader("MeshFormat")
        metaData = self.fileobj.readline().split()
        self.fileobj.seek(0)
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
     
    def _isolateData(self, title):
        """
        Gets all data between $[title] and $End[title], writes
        it out to its own file.
        """
        newF = tempfile.TemporaryFile()
        self._seekForHeader(title)
        
        # extract the actual data within section
        while True:
            line = self.fileobj.readline()
            if ("$End%s" % title) not in line: 
                newF.write(line) 
            else: break

        self.fileobj.seek(0) 
        newF.seek(0) # restore file positions
        return newF

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

        def formatForFiPy(arr): return arr.swapaxes(0,1)[::-1]

        allShapes  = nx.unique(shapeTypes).tolist()
        maxFaces   = max([self.numFacesPerCell[x] for x in allShapes])

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
            facesPerCell = self.numFacesPerCell[shapeType]

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

        return formatForFiPy(facesToVertices), formatForFiPy(cellsToFaces), facesDict
 
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

    def read(self):
        """
        0. Build cellsToVertices
        1. Recover needed vertexCoords and mapping from file using 
           cellsToVertices
        2. Build cellsToVertIDs proper from vertexCoords and vertex map
        3. Build faces
        4. Build cellsToFaces
        
        Isolate relevant data into three files, store in 
        `self.nodesFile` for $Nodes,
        `self.elemsFile` for $Elements. 
        `self.namesFile` for $PhysicalNames. 

        Returns vertexCoords, facesToVertexID, cellsToFaceID, 
                cellGlobalIDMap, ghostCellGlobalIDMap.
        """
        self.version, self.fileType, self.dataSize = self._getMetaData()
        self.nodesFile = self._isolateData("Nodes")
        self.elemsFile = self._isolateData("Elements")
        try:
            self.namesFile = self._isolateData("PhysicalNames")
        except EOFError, e:
            self.namesFile = None
            
        if self.dimensions is None:
            self.nodesFile.readline() # skip number of nodes

            # We assume we have a 2D file unless we find a node 
            # with a non-zero Z coordinate
            self.dimensions = 2
            for node in self.nodesFile:
                line   = node.split()
                
                newVert = [float(x) for x in line]
                
                if newVert[2] != 0.0:
                    self.dimensions = 3
                    break
                    
            self.nodesFile.seek(0)
            
        self.coordDimensions = self.coordDimensions or self.dimensions
            
        # we need a conditional here so we don't pick up 2D shapes in 3D
        if self.dimensions == 2: 
            self.numVertsPerFace = {1: 2} # line:       2 vertices
            self.numFacesPerCell = {2: 3, # triangle:   3 sides
                                    3: 4} # quadrangle: 4 sides
        elif self.dimensions == 3:
            self.numVertsPerFace = {2: 3, # triangle:   3 vertices
                                    3: 4} # quadrangle: 4 vertices
            self.numFacesPerCell = {4: 4, # tet:        4 sides
                                    5: 6} # hexahedron: 6 sides
        else:
            raise GmshException("Mesh has fewer than 2 or more than 3 dimensions")

        self.faceLenForShape = {2: 2, # triangle:   2 verts per face
                                3: 2, # quadrangle: 2 verts per face
                                4: 3, # tet:        3 verts per face
                                5: 4} # hexahedron: 4 verts per face

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
        cellsToVertIDs = self._translateVertIDToIdx(cellsToGmshVerts,
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
        faceEntitiesDict = dict()
        
        # translate Gmsh IDs to `vertexCoord` indices
        facesToVertIDs = self._translateVertIDToIdx(facesData.nodes,
                                                    vertIDtoIdx)
                                                    
        for face, physicalEntity, geometricalEntity in zip(facesToVertIDs, 
                                                           facesData.physicalEntities, 
                                                           facesData.geometricalEntities):
            faceEntitiesDict[' '.join([str(x) for x in sorted(face)])] = (physicalEntity, geometricalEntity)
            
        self.physicalFaceMap = nx.zeros(facesToV.shape[-1:])
        self.geometricalFaceMap = nx.zeros(facesToV.shape[-1:])
        for face in facesDict.keys():
            # not all faces are necessarily tagged
            if faceEntitiesDict.has_key(face):
                self.physicalFaceMap[facesDict[face]] = faceEntitiesDict[face][0]
                self.geometricalFaceMap[facesDict[face]] = faceEntitiesDict[face][1]
                
        self.physicalNames = self._parseNamesFile()
                                          
        parprint("Done with cells and faces.")
        return (vertexCoords, facesToV, cellsToF, 
                cellsData.idmap, ghostsData.idmap)
                
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
                if type(faceNum) not in [nx.int32, 
                                         nx.int64,
                                         nx.float32,
                                         nx.float64]:
                    continue
                for vertexNum in faceVertexIDs[..., faceNum]:
                    if vertexNum not in vertexList:
                        vertexList.append(vertexNum)
                        
            if dimensions == 2:
                vertexList = self._orderVertices(coords, vertexList)

            numVertices = len(vertexList)
            elementType = self._getElementType(numVertices, dimensions)
            self.fileobj.write("%s %s 0" % (str(i + 1), str(elementType)))

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
        for node in self.nodesFile:
            line   = node.split()
            nodeID = int(line[0])

            if nodeID == allVerts[nodeCount]:
                newVert = [float(x) for x in line[1:self.coordDimensions+1]]
                vertexCoords[nodeCount,:] = nx.array(newVert)
                nodeCount += 1

            if len(allVerts) == nodeCount: 
                break

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
        cellsData = _ElementData()
        ghostsData = _ElementData()
        facesData = _ElementData()
        
        cellOffset = -1 # this will be subtracted from gmsh ID to obtain global ID
        faceOffset = -1 # this will be subtracted from gmsh ID to obtain global ID
        pid = self.communicator.procID + 1
        
        self.elemsFile.readline() # skip number of elements
        for el in self.elemsFile:
            currLineInts = [int(x) for x in el.split()]
            elemType     = currLineInts[1]

            if elemType in self.numFacesPerCell.keys():
                # element is a cell
                if cellOffset == -1:
                    # if first valid shape
                    cellOffset = currLineInts[0]
                currLineInts[0] -= cellOffset

                numTags = currLineInts[2]
                tags    = currLineInts[3:(3+numTags)]
                
                physicalEntity = tags.pop(0)
                geometricalEntity = tags.pop(0)

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
                if faceOffset == -1:
                    # if first valid shape
                    faceOffset = currLineInts[0]
                currLineInts[0] -= faceOffset

                numTags = currLineInts[2]
                tags    = currLineInts[3:(3+numTags)]
                
                # the partition tags for faces don't seem to always be present 
                # and don't always make much sense when they are
                
                physicalEntity = tags.pop(0)
                geometricalEntity = tags.pop(0)

                facesData.add(currLine=currLineInts, elType=elemType, 
                              physicalEntity=physicalEntity, 
                              geometricalEntity=geometricalEntity)
                              
        self.elemsFile.close() # tempfile trashed

        return cellsData, ghostsData, facesData


    def _parseNamesFile(self):
        physicalNames = {
            1: dict(),
            2: dict(),
            3: dict()
        }
        if self.namesFile is not None:
            self.namesFile.readline() # skip number of elements
            for nm in self.namesFile:
                nm = nm.split()
                dim = int(nm.pop(0))
                num = int(nm.pop(0))
                name = " ".join(nm)[1:-1]
                physicalNames[dim][name] = int(num)
                
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
        >>> import subprocess
        >>> import tempfile
        
        >>> dir = tempfile.mkdtemp()
        
        >>> from fipy.meshes.grid2D import Grid2D
        >>> g = Grid2D(nx = 10, ny = 10)
        >>> f = openMSHFile(name=os.path.join(dir, "g.msh"), mode='w')
        >>> f.write(g)
        >>> f.close()

        >>> from fipy.meshes.uniformGrid2D import UniformGrid2D
        >>> ug = UniformGrid2D(nx = 10, ny = 10)
        >>> f = openMSHFile(name=os.path.join(dir, "ug.msh"), mode='w')
        >>> f.write(ug)
        >>> f.close()

        >>> from fipy.meshes import Tri2D
        >>> t = Tri2D(nx = 10, ny = 10)
        >>> f = openMSHFile(name=os.path.join(dir, "t.msh"), mode='w')
        >>> f.write(t)
        >>> f.close()

        >>> concat = ug + (t + ([10], [0]))
        >>> f = openMSHFile(name=os.path.join(dir, "concat.msh"), mode='w')
        >>> f.write(concat)
        >>> f.close()

        >>> from fipy.meshes import Grid3D
        >>> g3d = Grid3D(nx=10, ny=10, nz=30)
        >>> f = openMSHFile(name=os.path.join(dir, "g3d.msh"), mode='w')
        >>> f.write(g3d)
        >>> f.close()

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

def _makeMapVariables(mesh, physicalCellMap, geometricalCellMap, physicalFaceMap, 
                      geometricalFaceMap, physicalNames, faceDim, cellDim):
    """Utility function to make MeshVariables that define different domains in the mesh
    
    Separate function because Gmsh2D and Gmsh3D don't have a common ancestor and mixins are evil
    """
    from fipy.variables.cellVariable import CellVariable
    from fipy.variables.faceVariable import FaceVariable

    physicalCellMap = CellVariable(mesh=mesh, value=physicalCellMap)
    geometricalCellMap = CellVariable(mesh=mesh, value=geometricalCellMap)
    physicalFaceMap = FaceVariable(mesh=mesh, value=physicalFaceMap)
    geometricalFaceMap = FaceVariable(mesh=mesh, value=geometricalFaceMap)

    physicalCells = dict()
    for name in physicalNames[cellDim].keys():
        physicalCells[name] = physicalCellMap == physicalNames[cellDim][name]
        
    physicalFaces = dict()
    for name in physicalNames[faceDim].keys():
        physicalFaces[name] = physicalFaceMap == physicalNames[faceDim][name]
        
    return (physicalCellMap, geometricalCellMap, physicalCells,
            physicalFaceMap, geometricalFaceMap, physicalFaces)

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
    ... ''' % locals())

    It can be easier to specify certain domains and boundaries within Gmsh 
    than it is to define the same domains and boundaries with FiPy expressions.
    
    Here we compare obtaining the same Cells and Faces using FiPy's 
    parametric descriptions and Gmsh's labels.
    
    >>> x, y = squaredCircle.cellCenters

    >>> middle = ((x**2 + y**2 <= radius**2) 
    ...           & ~((x > -side/2) & (x < side/2)
    ...               & (y > -side/2) & (y < side/2)))

    >>> print (middle == squaredCircle.physicalCells["Middle"]).all()
    True
    
    >>> X, Y = squaredCircle.faceCenters

    >>> NW = ((X**2 + Y**2 > (1.99*radius)**2) 
    ...       & (X**2 + Y**2 < (2.01*radius)**2)
    ...       & (X <= 0) & (Y >= 0))

    >>> print (NW == squaredCircle.physicalFaces["NW"]).all()
    True
    """
    
    def __init__(self, 
                 arg, 
                 coordDimensions=2, 
                 communicator=parallel, 
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
         self.gCellGlobalIDs) = self.mshFile.read()
         
        self.mshFile.close()

        if communicator.Nproc > 1:
            self.globalNumberOfCells = communicator.sumAll(len(self.cellGlobalIDs))
            parprint("  I'm solving with %d cells total." % self.globalNumberOfCells)
            parprint("  Got global number of cells")

        Mesh2D.__init__(self, vertexCoords=verts,
                              faceVertexIDs=faces,
                              cellFaceIDs=cells,
                              communicator=communicator)
                   
        (self.physicalCellMap,
         self.geometricalCellMap,
         self.physicalCells,
         self.physicalFaceMap,
         self.geometricalFaceMap,
         self.physicalFaces) = self.mshFile.makeMapVariables(mesh=self)
        
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
        ... ''')

        >>> print circ.cellVolumes[0] > 0
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
        self.mshFile  = openMSHFile(arg, 
                                    dimensions=3, 
                                    communicator=communicator,
                                    order=order,
                                    mode='r')

        (verts,
         faces,
         cells,
         self.cellGlobalIDs,
         self.gCellGlobalIDs) = self.mshFile.read()
         
        self.mshFile.close()

        Mesh.__init__(self, vertexCoords=verts,
                            faceVertexIDs=faces,
                            cellFaceIDs=cells,
                            communicator=communicator)

        if self.communicator.Nproc > 1:
            self.globalNumberOfCells = self.communicator.sumAll(len(self.cellGlobalIDs))
   
        (self.physicalCellMap,
         self.geometricalCellMap,
         self.physicalCells,
         self.physicalFaceMap,
         self.geometricalFaceMap,
         self.physicalFaces) = self.mshFile.makeMapVariables(mesh=self)
        
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
                 coordDimensions=2, communicator=parallel, order=1,
                 background=None):
        self.dx = dx
        self.dy = dy or dx
        self.nx = nx
        self.ny = ny or nx

        arg = self._makeGridGeo(self.dx, self.dy, self.nx, self.ny)

        Gmsh2D.__init__(self, arg, coordDimensions, communicator, order, background)
    
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
    
    from fipy.meshes import Grid2D
    from fipy.meshes import Tri2D
    from fipy.meshes import Grid3D
    from fipy.meshes import CylindricalGrid2D
    from fipy.meshes import Gmsh2D
    from fipy.meshes import GmshGrid2D
    from fipy.variables.cellVariable import CellVariable
    
    import tempfile
    import subprocess
    
    dir = tempfile.mkdtemp()

    a = Grid2D(dx = 1.0, dy = 1.0, nx = 10, ny = 10)
    a2 = Grid2D(dx = 1.0, dy = 1.0, nx = 10, ny = 10) + ([10], [0])
    
    x, y = a.cellCenters
    avar = CellVariable(mesh=a, name="binky", value=1. / (x+y))
    
    f1 = openMSHFile(name=os.path.join(dir, "a.msh"), mode='w')
    f1.write(a)
    f1.write(avar)
    f1.close()
    
    subprocess.Popen(["gmsh", os.path.join(dir, "a.msh")])
    raw_input("Grid2D... Press enter.")
    
    a_ref = GmshGrid2D(dx=1., dy=1., nx=10, ny=10, background=avar)
    
    f1 = openMSHFile(name=os.path.join(dir, "a_ref.msh"), mode='w')
    f1.write(a_ref)
    f1.close()

    subprocess.Popen(["gmsh", os.path.join(dir, "a_ref.msh")])
    raw_input("Refined Grid2D... Press enter.")

    b = Tri2D(dx = 1.0, dy = 1.0, nx = 10, ny = 10)
    
    f2 = openMSHFile(name=os.path.join(dir, "b.msh"), mode='w')
    f2.write(b)
    f2.close()

    c = Grid3D(dx = 1.0, dy = 1.0, nx = 20, ny = 20, nz = 40)
    
    f3 = openMSHFile(name=os.path.join(dir, "c.msh"), mode='w')
    f3.write(c)
    f3.close()

    subprocess.Popen(["gmsh", os.path.join(dir, "c.msh")])
    raw_input("Grid3D... Press enter.")

    d = a + a2
    f4 = openMSHFile(name=os.path.join(dir, "d.msh"), mode='w')
    f4.write(d)
    f4.close()

    subprocess.Popen(["gmsh", os.path.join(dir, "d.msh")])
    raw_input("Concatenated grid... Press enter.")
 
    e = a + (b + ([0], [10]))
    f5 = openMSHFile(name=os.path.join(dir, "e.msh"), mode='w')
    f5.write(e)
    f5.close()

    subprocess.Popen(["gmsh", os.path.join(dir, "e.msh")])
    raw_input("Tri2D + Grid2D... Press enter.")

    cyl = CylindricalGrid2D(nx = 10, ny = 10)
    f6 = openMSHFile(name=os.path.join(dir, "cyl.msh"), mode='w')
    f6.write(cyl)
    f6.close()
    
    subprocess.Popen(["gmsh", os.path.join(dir, "cyl.msh")])
    raw_input("CylindricalGrid2D... Press enter.")

    circle = Gmsh2D('''
         cellSize = 0.05;
         radius = 1;
         Point(1) = {0, 0, 0, cellSize};
         Point(2) = {-radius, 0, 0, cellSize};
         Point(3) = {0, radius, 0, cellSize};
         Point(4) = {radius, 0, 0, cellSize};
         Point(5) = {0, -radius, 0, cellSize};
         Circle(6) = {2, 1, 3};
         Circle(7) = {3, 1, 4};
         Circle(8) = {4, 1, 5};
         Circle(9) = {5, 1, 2};
         Line Loop(10) = {6, 7, 8, 9};
         Plane Surface(11) = {10};
         Recombine Surface{11};
    ''')   
    f7 = openMSHFile(name=os.path.join(dir, "cir.msh"), mode='w')
    f7.write(circle)
    f7.close()
    
    subprocess.Popen(["gmsh", os.path.join(dir, "cir.msh")])
    raw_input("Circle... Press enter.")
    
    import shutil
    shutil.rmtree(dir)
