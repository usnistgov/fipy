#!/usr/bin/env python

##
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "mesh1D.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

"""Generic mesh class using numerix to do the calculations

    Meshes contain cells, faces, and vertices.

    This is built for a non-mixed element mesh.
"""

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools import serialComm

from fipy.meshes.mesh import Mesh
from fipy.meshes.representations.meshRepresentation import _MeshRepresentation
from fipy.meshes.topologies.meshTopology import _Mesh1DTopology

__all__ = ["Mesh1D"]

class Mesh1D(Mesh):
    def __init__(self, vertexCoords, faceVertexIDs, cellFaceIDs, communicator=serialComm, _RepresentationClass=_MeshRepresentation, _TopologyClass=_Mesh1DTopology):
        super(Mesh1D, self).__init__(vertexCoords=vertexCoords, faceVertexIDs=faceVertexIDs, cellFaceIDs=cellFaceIDs, communicator=communicator,
                                     _RepresentationClass=_RepresentationClass, _TopologyClass=_TopologyClass)

    def _calcScaleArea(self):
        return 1.

    def _calcScaleVolume(self):
        return self.scale['length']

    def _calcFaceAreas(self):
        return numerix.ones(self.numberOfFaces, 'd')

    def _calcFaceNormals(self):
        faceNormals = numerix.array((numerix.ones(self.numberOfFaces, 'd'),))
        # The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
        if self.numberOfFaces > 0:
            faceNormals[...,0] = -faceNormals[...,0]
        return faceNormals

    def _calcFaceTangents(self):
        faceTangents1 = numerix.zeros(self.numberOfFaces, 'd')[numerix.NewAxis, ...]
        faceTangents2 = numerix.zeros(self.numberOfFaces, 'd')[numerix.NewAxis, ...]
        return faceTangents1, faceTangents2

    def _translate(self, vector):
        newCoords = self.vertexCoords + vector
        newmesh = Mesh1D(newCoords, numerix.array(self.faceVertexIDs), numerix.array(self.cellFaceIDs))
        return newmesh

    def __mul__(self, factor):
        newCoords = self.vertexCoords * factor
        newmesh = Mesh1D(newCoords, numerix.array(self.faceVertexIDs), numerix.array(self.cellFaceIDs))
        return newmesh

    @property
    def _VTKCellType(self):
        try:
            from tvtk.api import tvtk
        except ImportError, e:
            from enthought.tvtk.api import tvtk
        return tvtk.Line().cell_type
