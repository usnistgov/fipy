#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh2D.py"
 #                                    created: 11/10/03 {2:44:42 PM} 
 #                                last update: 6/2/04 {3:56:32 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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
 ##

"""Generic mesh class using Numeric to do the calculations

    Meshes contain cells, faces, and vertices.

    This is built for a non-mixed element mesh.
"""

import Numeric
import MA
from fipy.meshes.numMesh.mesh import Mesh
import fipy.tools.array
import fipy.tools.vector as vector


class Mesh2D(Mesh):
    def calcFaceAreas(self):
        faceVertexCoords = Numeric.take(self.vertexCoords, self.faceVertexIDs)
        tangent = faceVertexCoords[:,1] - faceVertexCoords[:,0]
        self.faceAreas = fipy.tools.array.sqrtDot(tangent, tangent)

    def calcFaceNormals(self):
        faceVertexCoords = Numeric.take(self.vertexCoords, self.faceVertexIDs)
        t1 = faceVertexCoords[:,1,:] - faceVertexCoords[:,0,:]
        self.faceNormals = t1.copy()
        self.faceNormals[:,0] = -t1[:,1] / Numeric.sqrt(t1[:,1]**2 + t1[:,0]**2)
        self.faceNormals[:,1] = t1[:,0] / Numeric.sqrt(t1[:,1]**2 + t1[:,0]**2)

    def calcFaceTangents(self):
        tmp = Numeric.transpose(Numeric.array((-self.faceNormals[:,1], self.faceNormals[:,0])))
        mag = fipy.tools.array.sqrtDot(tmp, tmp)
        self.faceTangents1 = tmp / mag[:,Numeric.NewAxis]
        self.faceTangents2 = Numeric.zeros(self.faceTangents1.shape, 'd')

    def calcHigherOrderScalings(self):
	self.scale['area'] = self.scale['length']
	self.scale['volume'] = self.scale['length']**2

    def translate(self, vector):
        newCoords = self.vertexCoords + vector
        newmesh = Mesh2D(newCoords, Numeric.array(self.faceVertexIDs), Numeric.array(self.cellFaceIDs))
        return newmesh

    def dilate(self, factor):
        newCoords = self.vertexCoords * factor
        newmesh = Mesh2D(newCoords, Numeric.array(self.faceVertexIDs), Numeric.array(self.cellFaceIDs))
        return newmesh

    def meshAdd(self, other):
        a = self.getAddedMeshValues(other)
        return Mesh2D(a[0], a[1], a[2])

    def getNonOrthogonality(self):
        exteriorFaceArray = Numeric.zeros((self.faceCellIDs.shape[0],))
        Numeric.put(exteriorFaceArray, self.getExteriorFaceIDs(), 1)
        unmaskedFaceCellIDs = MA.filled(self.faceCellIDs, 0) ## what we put in for the "fill" doesn't matter because only exterior faces have anything masked, and exterior faces have their displacement vectors set to zero.
        ## if it's an exterior face, make the "displacement vector" equal to zero so the cross product will be zero.
        faceDisplacementVectors = Numeric.where(Numeric.array(zip(exteriorFaceArray, exteriorFaceArray)), 0.0, Numeric.take(self.getCellCenters(), unmaskedFaceCellIDs[:, 1]) - Numeric.take(self.getCellCenters(), unmaskedFaceCellIDs[:, 0]))
        faceCrossProducts = (faceDisplacementVectors[:, 0] * self.faceNormals[:, 1]) - (faceDisplacementVectors[:, 1] * self.faceNormals[:, 0])
        faceDisplacementVectorLengths = Numeric.maximum(((faceDisplacementVectors[:, 0] ** 2) + (faceDisplacementVectors[:, 1] ** 2)) ** 0.5, 1.e-100)
        faceWeightedNonOrthogonalities = abs(faceCrossProducts / faceDisplacementVectorLengths) * self.faceAreas
        cellFaceWeightedNonOrthogonalities = Numeric.take(faceWeightedNonOrthogonalities, self.cellFaceIDs)
        cellFaceAreas = Numeric.take(self.faceAreas, self.cellFaceIDs)
        cellTotalWeightedValues = Numeric.add.reduce(cellFaceWeightedNonOrthogonalities, axis = 1)
        cellTotalFaceAreas = Numeric.add.reduce(cellFaceAreas, axis = 1)
        return (cellTotalWeightedValues / cellTotalFaceAreas)
