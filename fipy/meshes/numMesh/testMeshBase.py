#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "testMeshBase.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 3/5/04 {1:54:20 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
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
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

"""Test numeric implementation of the mesh
"""
 
import unittest
from fivol.tests.testBase import TestBase
import fivol.tests.testProgram
import Numeric
import MA

from mesh import Mesh

class TestMeshBase(TestBase):
    def testResult(self):
        pass

    def testExteriorFaces(self):
        self.assertEqual(self.externalFaces, [face.getID() for face in self.mesh.getExteriorFaces()])
        
    def testInternalFaces(self):
        self.assertArrayEqual(self.internalFaces, [face.getID() for face in self.mesh.getInteriorFaces()])

    def testFaceCellIds(self):
        self.assertArrayEqual(self.faceCellIds, self.mesh.getFaceCellIDs())

    def testFaceAreas(self):
        self.assertArrayWithinTolerance(self.faceAreas, self.mesh.getFaceAreas())
        
    def testFaceCenters(self):
        self.assertArrayWithinTolerance(self.faceCenters, self.mesh.getFaceCenters())

    def testFaceNormals(self):
        self.assertArrayWithinTolerance(self.faceNormals, self.mesh.getFaceNormals())

    def testCellToFaceOrientations(self):
	self.assertArrayEqual(self.cellToFaceOrientations, self.mesh.getCellFaceOrientations())

    def testCellVolumes(self):
	self.assertArrayWithinTolerance(self.cellVolumes, self.mesh.getCellVolumes())

    def testCellCenters(self):
	self.assertArrayWithinTolerance(self.cellCenters, self.mesh.getCellCenters())

    def testFaceToCellDistances(self):
	self.assertArrayWithinTolerance(self.faceToCellDistances, self.mesh.getFaceToCellDistances())

    def testCellDistances(self):
	self.assertArrayWithinTolerance(self.cellDistances, self.mesh.getCellDistances())

    def testFaceToCellDistanceRatios(self):
        self.assertArrayWithinTolerance(self.faceToCellDistanceRatios, self.mesh.getFaceToCellDistanceRatio())

    def testAreaProjections(self):
        self.assertArrayWithinTolerance(self.areaProjections, self.mesh.getAreaProjections())

    def testTangents1(self):
        self.assertArrayWithinTolerance(self.tangents1, self.mesh.getFaceTangents1())

    def testTangents2(self):
        self.assertArrayWithinTolerance(self.tangents2, self.mesh.getFaceTangents2())
