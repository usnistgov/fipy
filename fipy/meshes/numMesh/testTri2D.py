#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:00:01 PM} 
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
from fipy.tests.testBase import TestBase
import fipy.tests.testProgram
import Numeric
from fipy.meshes.numMesh.tri2D import Tri2D
import MA
from fipy.meshes.numMesh.testMeshBase import TestMeshBase
import fipy.tools.dump as dump

class TestTri2D(TestMeshBase):
    def setUp(self):
	dx = 0.5
	dy = 2.
	nx = 3
	ny = 2
	
        self.mesh = Tri2D(nx = nx, ny = ny, dx = dx, dy = dy)     
        
        self.vertices = Numeric.array(((0.0, 0.0), (0.5, 0.0), (1.0, 0.0), (1.5, 0.0),
                                       (0.0, 2.0), (0.5, 2.0), (1.0, 2.0), (1.5, 2.0),
                                       (0.0, 4.0), (0.5, 4.0), (1.0, 4.0), (1.5, 4.0),
                                       (0.25, 1.0), (0.75, 1.0), (1.25, 1.0),
                                       (0.25, 3.0), (0.75, 3.0), (1.25, 3.0)))
        self.faces = Numeric.array(((1, 0), (2, 1), (3, 2),
                                    (4, 5), (5, 6), (6, 7),
                                    (8, 9), (9, 10), (10, 11),
                                    (0, 4), (5, 1), (6, 2), (7, 3),
                                    (4, 8), (9, 5), (10, 6), (11, 7),
                                    (12, 0), (13, 1), (14, 2), (15, 4), (16, 5), (17, 6),
                                    (1, 12), (2, 13), (3, 14), (5, 15), (6, 16), (7, 17),
                                    (12, 4), (13, 5), (14, 6), (15, 8), (16, 9), (17, 10),
                                    (12, 5), (13, 6), (14, 7), (15, 9), (16, 10), (17, 11)))
        self.cells = Numeric.array(((10, 35, 23), (11, 36, 24), (12, 37, 25), (14, 38, 26), (15, 39, 27), (16, 40, 28),
                                    (3, 29, 35), (4, 30, 36), (5, 31, 37), (6, 32, 38), (7, 33, 39), (8, 34, 40),
                                    (9, 17, 29), (10, 18, 30), (11, 19, 31), (13, 20, 32), (14, 21, 33), (15, 22, 34),
                                    (0, 23, 17), (1, 24, 18), (2, 25, 19), (3, 26, 20), (4, 27, 21), (5, 28, 22)))


        self.externalFaces = Numeric.array((0, 1, 2, 6, 7, 8, 9 , 12, 13, 16))
        self.internalFaces = Numeric.array((3, 4, 5, 10, 11, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40))

        
        self.faceCellIds = MA.masked_values(((18, -1), (19, -1), (20, -1), (6, 21), (7, 22), (8, 23), (9, -1), (10, -1), (11, -1),
                                             (12, -1), (0, 13), (1, 14), (2, -1), (15, -1), (3, 16), (4, 17), (5, -1),
                                             (12, 18), (13, 19), (14, 20), (15, 21), (16, 22), (17, 23),
                                             (0, 18), (1, 19), (2, 20), (3, 21), (4, 22), (5, 23),
                                             (6, 12), (7, 13), (8, 14), (9, 15), (10, 16), (11, 17),
                                             (0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)), -1)
        
        d = (Numeric.sqrt((dx*dx)+(dy*dy))) / 2.0 ## length of diagonal edges  

        self.faceAreas = Numeric.array((0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                                        2, 2, 2, 2, 2, 2, 2, 2,
                                        d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d))
                                        
        
        faceCoords = Numeric.take(self.vertices, self.faces)
        self.faceCenters = (faceCoords[:,0] + faceCoords[:,1]) / 2.

        xc = dy  / Numeric.sqrt((dx * dx) + (dy * dy))
        yc = dx  / Numeric.sqrt((dx * dx) + (dy * dy))
        
        self.faceNormals = Numeric.array(((0.0, -1.0), (0.0, -1.0), (0.0, -1.0),
                                          (0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
                                          (0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
                                          (-1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0),
                                          (-1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0),
                                          (xc, -yc), (xc, -yc), (xc, -yc), (xc, -yc), (xc, -yc), (xc, -yc),
                                          (-xc, -yc), (-xc, -yc), (-xc, -yc), (-xc, -yc), (-xc, -yc), (-xc, -yc),
                                          (-xc, -yc), (-xc, -yc), (-xc, -yc), (-xc, -yc), (-xc, -yc), (-xc, -yc),
                                          (-xc, yc), (-xc, yc), (-xc, yc), (-xc, yc), (-xc, yc), (-xc, yc)))
        
	self.cellToFaceOrientations = Numeric.array(((1, 1, 1), (1, 1, 1), (1, 1, 1), (1, 1, 1), (1, 1, 1), (1, 1, 1),
                                                     (1, 1, -1), (1, 1, -1), (1, 1, -1), (1, 1, -1), (1, 1, -1), (1, 1, -1),
                                                     (1, 1, -1), (-1, 1, -1), (-1, 1, -1), (1, 1, -1), (-1, 1, -1), (-1, 1, -1),
                                                     (1, -1, -1), (1, -1, -1), (1, -1, -1), (-1, -1, -1), (-1, -1, -1), (-1, -1, -1)))
					 
        self.cellVolumes = Numeric.array((0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                                       0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                                       0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25))

        sixth = 1.0 / 6.0
        
	self.cellCenters = Numeric.array(((5*sixth, 0.5), (11*sixth, 0.5), (17*sixth, 0.5), (5*sixth, 1.5), (11*sixth, 1.5), (17*sixth, 1.5),
                                         (0.5, 5*sixth), (1.5, 5*sixth), (2.5, 5*sixth), (0.5, 11*sixth), (1.5, 11*sixth), (2.5, 11*sixth),
                                         (1*sixth, 0.5), (7*sixth, 0.5), (13*sixth, 0.5), (1*sixth, 1.5), (7*sixth, 1.5), (13*sixth, 1.5),
                                         (0.5, 1*sixth), (1.5, 1*sixth), (2.5, 1*sixth), (0.5, 7*sixth), (1.5, 7*sixth), (2.5, 7*sixth)))

        self.cellCenters = self.cellCenters * Numeric.array((dx, dy))


        yd = Numeric.sqrt( ((dx/12.0)*(dx/12.0)) + ((dy/4.0)*(dy/4.0)) )
        xd = Numeric.sqrt( ((dx/4.0)*(dx/4.0)) + ((dy/12.0)*(dy/12.0)) )
        
	self.faceToCellDistances = MA.masked_values(((dy/6.0, -1), (dy/6.0, -1), (dy/6.0, -1), (dy/6.0, dy/6.0), (dy/6.0, dy/6.0), (dy/6.0, dy/6.0), (dy/6.0, -1), (dy/6.0, -1), (dy/6.0, -1),
                                                     (dx/6.0, -1), (dx/6.0, dx/6.0), (dx/6.0, dx/6.0), (dx/6.0, -1), (dx/6.0, -1), (dx/6.0, dx/6.0), (dx/6.0, dx/6.0), (dx/6.0, -1),
                                                     (yd, xd), (yd, xd), (yd, xd), (yd, xd), (yd, xd), (yd, xd),
                                                     (yd, xd), (yd, xd), (yd, xd), (yd, xd), (yd, xd), (yd, xd),
                                                     (xd, yd), (xd, yd), (xd, yd), (xd, yd), (xd, yd), (xd, yd),
                                                     (yd, xd), (yd, xd), (yd, xd), (yd, xd), (yd, xd), (yd, xd)), -1)

        dd = Numeric.sqrt((dx*dx)+(dy*dy)) / 3.0
        
	self.cellDistances = Numeric.array((dy/6.0, dy/6.0, dy/6.0, dy/3.0, dy/3.0, dy/3.0, dy/6.0, dy/6.0, dy/6.0,
                                            dx/6.0, dx/3.0, dx/3.0, dx/6.0, dx/6.0, dx/3.0, dx/3.0, dx/6.0,
                                            dd, dd, dd, dd, dd, dd, dd, dd, dd, dd, dd, dd,
                                            dd, dd, dd, dd, dd, dd, dd, dd, dd, dd, dd, dd))
            
        self.faceToCellDistanceRatios = self.faceToCellDistances[:,0] / self.cellDistances

        self.areaProjections = self.faceNormals * self.faceAreas[:,Numeric.NewAxis]

        self.facesBottom = (0, 1, 2)
	self.facesTop = (6, 7, 8)
	self.facesLeft = (9, 13)
	self.facesRight = (12, 16)


        self.tangents1 = Numeric.array(((1.0, 0.0), (1.0, 0.0), (1.0, 0.0),
                                        (-1.0, 0.0), (-1.0, 0.0), (-1.0, 0.0),
                                        (-1.0, 0.0), (-1.0, 0.0), (-1.0, 0.0),
                                        (0.0, -1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
                                        (0.0, -1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
                                        (yc, xc),(yc, xc),(yc, xc),(yc, xc),(yc, xc),(yc, xc),
                                        (yc, -xc),(yc, -xc),(yc, -xc),(yc, -xc),(yc, -xc),(yc, -xc),
                                        (yc, -xc),(yc, -xc),(yc, -xc),(yc, -xc),(yc, -xc),(yc, -xc),
                                        (-yc, -xc),(-yc, -xc),(-yc, -xc),(-yc, -xc),(-yc, -xc),(-yc, -xc)))
                                        

        self.tangents2 = Numeric.zeros((41, 2))

        self.cellToCellIDs = MA.masked_values(((13, 18, 6), (14, 19, 7), (-1, 20, 8), (16, 21, 9), (17, 22, 10), (-1, 23, 11),
                                               (21, 12, 0), (22, 13, 1), (23, 14, 2), (-1, 15, 3), (-1, 16, 4), (-1, 17, 5),
                                               (-1, 18, 6), (0, 19, 7), (1, 20, 8), (-1, 21, 9), (3, 22, 10), (4, 23, 11),
                                               (-1, 12, 0), (-1, 13, 1), (-1, 14, 2), (6, 15, 3), (7, 16, 4), (8, 17, 5)), -1)


        ## self.cellToCellDistances = Numeric.take(self.cellDistances, self.cells)
        self.cellToCellDistances = Numeric.array([[dx/3.0, dd, dd],
                                    [ dx/3.0, dd, dd],
                                    [ dx/6.0, dd, dd],
                                    [ dx/3.0, dd, dd],
                                    [ dx/3.0, dd, dd],
                                    [ dx/6.0, dd, dd],
                                    [ dy/3.0, dd, dd],
                                    [ dy/3.0, dd, dd],
                                    [ dy/3.0, dd, dd],
                                    [ dy/6.0, dd, dd],
                                    [ dy/6.0, dd, dd],
                                    [ dy/6.0, dd, dd],
                                    [ dx/6.0, dd, dd],
                                    [ dx/3.0, dd, dd],
                                    [ dx/3.0, dd, dd],
                                    [ dx/6.0, dd, dd],
                                    [ dx/3.0, dd, dd],
                                    [ dx/3.0, dd, dd],
                                    [ dy/6.0, dd, dd],
                                    [ dy/6.0, dd, dd],
                                    [ dy/6.0, dd, dd],
                                    [ dy/3.0, dd, dd],
                                    [ dy/3.0, dd, dd],
                                    [ dy/3.0, dd, dd]])

        self.interiorCellIDs = Numeric.array((0, 1, 3, 4, 6, 7, 8, 13, 14, 16, 17, 21, 22, 23))

        self.exteriorCellIDs = Numeric.array((2, 5, 9, 10, 11, 12, 15, 18, 19, 20))

                                       
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

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestTri2D))
    return theSuite

if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')


