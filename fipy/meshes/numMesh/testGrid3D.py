#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 9/3/04 {10:41:44 PM} 
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
from fipy.tests.testBase import TestBase
import fipy.tests.testProgram
import Numeric
from fipy.meshes.numMesh.grid3D import Grid3D
import MA
from fipy.meshes.numMesh.testMeshBase import TestMeshBase
import fipy.tools.dump as dump

class TestGrid(TestMeshBase):
    def setUp(self):
	dx = 0.5
	dy = 2.
        dz = 4.
	nx = 3
	ny = 2
        nz = 1
	
        self.mesh = Grid3D(nx = nx, ny = ny, nz = nz, dx = dx, dy = dy, dz = dz)     
        
        self.vertices = Numeric.array(((0., 0., 0.), (1., 0., 0.), (2., 0., 0.), (3., 0., 0.),
                                       (0., 1., 0.), (1., 1., 0.), (2., 1., 0.), (3., 1., 0.),
                                       (0., 2., 0.), (1., 2., 0.), (2., 2., 0.), (3., 2., 0.),
                                       (0., 0., 1.), (1., 0., 1.), (2., 0., 1.), (3., 0., 1.),
                                       (0., 1., 1.), (1., 1., 1.), (2., 1., 1.), (3., 1., 1.),
                                       (0., 2., 1.), (1., 2., 1.), (2., 2., 1.), (3., 2., 1.)))
        
        self.vertices = self.vertices * Numeric.array((dx, dy, dz))
        
        self.faces = Numeric.array(((0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (4, 5, 9, 8), (5, 6, 10, 9), (6, 7, 11, 10),
                                    (12, 13, 17, 16), (13, 14, 18, 17), (14, 15, 19, 18), (16, 17, 21, 20), (17, 18, 22, 21), (18, 19, 23, 22),
                                    (0, 1, 13, 12), (1, 2, 14, 13), (2, 3, 15, 14), (4, 5, 17, 16), (5, 6, 18, 17), (6, 7, 19, 18), (8, 9, 21, 20), (9, 10, 22, 21), (10, 11, 23, 22),
                                    (0, 4, 16, 12), (1, 5, 17, 13), (2, 6, 18, 14), (3, 7, 19, 15), (4, 8, 20, 16), (5, 9, 21, 17), (6, 10, 22, 18), (7, 11, 23, 19)))
            
        self.cells = Numeric.array(((21, 22, 12, 15, 0, 6),
                                   (22, 23, 13, 16, 1, 7),
                                   (23, 24, 14, 17, 2, 8),
                                   (25, 26, 15, 18, 3, 9),
                                   (26, 27, 16, 19, 4, 10),
                                   (27, 28, 17, 20, 5, 11)))


        self.externalFaces = Numeric.array((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 24, 25, 28))
        self.internalFaces = Numeric.array((15, 16, 17, 22, 23, 26, 27))
        
        self.faceCellIds = MA.masked_values(((0, -1), (1, -1), (2, -1), (3, -1), (4, -1), (5, -1),
                                             (0, -1), (1, -1), (2, -1), (3, -1), (4, -1), (5, -1),
                                             (0, -1), (1, -1), (2, -1), (0, 3), (1, 4), (2, 5), (3, -1), (4, -1), (5, -1),
                                             (0, -1), (0, 1), (1, 2), (2, -1), (3, -1), (3, 4), (4, 5), (5, -1)), -1)
        xy = dx * dy
        xz = dx * dz
        yz = dy * dz
        
        self.faceAreas = Numeric.array((xy, xy, xy, xy, xy, xy, xy, xy, xy, xy, xy, xy,
                                        xz, xz, xz, xz, xz, xz, xz, xz, xz,
                                        yz, yz, yz, yz, yz, yz, yz, yz))
        
        faceCoords = Numeric.take(self.vertices, self.faces)

        self.faceCenters = (faceCoords[:,0] + faceCoords[:,1] + faceCoords[:,2] + faceCoords[:, 3]) / 4.

        self.faceNormals = Numeric.array(((0, 0, -1), (0, 0, -1), (0, 0, -1), (0, 0, -1), (0, 0, -1), (0, 0, -1),
                                          (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1),
                                          (0, -1, 0), (0, -1, 0), (0, -1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0),
                                          (-1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (-1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0)))
        
	self.cellToFaceOrientations = Numeric.array(((1, 1, 1, 1, 1, 1),
                                                    (-1, 1, 1, 1, 1, 1),
                                                    (-1, 1, 1, 1, 1, 1),
                                                    (1, 1, -1, 1, 1, 1),
                                                    (-1, 1, -1, 1, 1, 1),
                                                    (-1, 1, -1, 1, 1, 1)))
					 
	self.cellVolumes = Numeric.array((dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz, dx*dy*dz))

	self.cellCenters = Numeric.array(((dx/2.,dy/2.,dz/2.), (3.*dx/2.,dy/2.,dz/2.), (5.*dx/2.,dy/2.,dz/2.),
                                          (dx/2.,3.*dy/2.,dz/2.), (3.*dx/2.,3.*dy/2.,dz/2.), (5.*dx/2.,3.*dy/2.,dz/2.)))
					  
	self.faceToCellDistances = MA.masked_values(((dz/2, -1), (dz/2, -1), (dz/2, -1), (dz/2, -1), (dz/2, -1), (dz/2, -1),
                                                     (dz/2, -1), (dz/2, -1), (dz/2, -1), (dz/2, -1), (dz/2, -1), (dz/2, -1),
                                                     (dy/2, -1), (dy/2, -1), (dy/2, -1), (dy/2, dy/2), (dy/2, dy/2), (dy/2, dy/2), (dy/2, -1), (dy/2, -1), (dy/2, -1),
                                                     (dx/2, -1), (dx/2, dx/2), (dx/2, dx/2), (dx/2, -1), (dx/2, -1), (dx/2, dx/2), (dx/2, dx/2), (dx/2, -1)), -1) 
                                                     
					  
	self.cellDistances = Numeric.array((dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2, dz/2,
                                            dy/2, dy/2, dy/2, dy, dy, dy, dy/2, dy/2, dy/2,
                                            dx/2, dx, dx, dx/2, dx/2, dx, dx, dx/2))
        
        self.faceToCellDistanceRatios = self.faceToCellDistances[:,0] / self.cellDistances

        self.areaProjections = self.faceNormals * self.faceAreas[:,Numeric.NewAxis]

        self.facesFront = (0, 1, 2, 3, 4, 5)
	self.facesBack = (6, 7, 8, 9, 10, 11)
	self.facesBottom = (12, 13, 14)
	self.facesTop = (18, 19, 20)
        self.facesLeft = (21, 25)
        self.facesRight = (24, 28)

        self.tangents1 = Numeric.array(((1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0),
                                        (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0),
                                        (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0)))

        self.tangents2 = Numeric.array(((0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0),
                                        (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1),
                                        (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1), (0, 0, 1)))

        self.cellToCellIDs = MA.masked_values(((-1, -1, -1, 3, -1, 1),
                                               (-1, -1, -1,  4, 0, 2),
                                               (-1, -1, -1, 5, 1, -1),
                                               (-1, -1, 0, -1, -1, 4),
                                               (-1, -1, 1, -1, 3, 5),
                                               (-1 , -1, 2, -1, 4, -1)), -1)

        self.cellToCellDistances = Numeric.take(self.cellDistances, self.cells)
                                                   
        self.interiorCellIDs = Numeric.array(())

        self.exteriorCellIDs = Numeric.array((0, 1, 2, 3, 4, 5))

                                       
    def testVertices(self):
        self.assertArrayEqual(self.vertices, self.mesh.createVertices())

    def testFaces(self):
        self.assertArrayEqual(self.faces, self.mesh.createFaces())

    def testCells(self):
	self.assertArrayEqual(self.cells, self.mesh.createCells())
	
    def testFacesBottom(self):
	self.assertArrayEqual(self.facesBottom, [face.getID() for face in self.mesh.getFacesBottom()])
	
    def testFacesTop(self):
	self.assertArrayEqual(self.facesTop, [face.getID() for face in self.mesh.getFacesTop()])

    def testFacesLeft(self):
	self.assertArrayEqual(self.facesLeft, [face.getID() for face in self.mesh.getFacesLeft()])

    def testFacesRight(self):
	self.assertArrayEqual(self.facesRight, [face.getID() for face in self.mesh.getFacesRight()])

    def testFacesFront(self):
        self.assertArrayEqual(self.facesFront, [face.getID() for face in self.mesh.getFacesFront()])

    def testFacesBack(self):
        self.assertArrayEqual(self.facesBack, [face.getID() for face in self.mesh.getFacesBack()])

class TestGridPickle(TestGrid):
    def setUp(self):
        TestGrid.setUp(self)
        import tempfile
        import os
        tmp = tempfile.gettempdir()
        fileName = os.path.join(tmp, 'data')
        pickledMesh = dump.write(self.mesh, fileName)
        self.mesh = dump.read(fileName)

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestGrid))
    theSuite.addTest(unittest.makeSuite(TestGridPickle))
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')
