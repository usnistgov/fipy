#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 11/16/04 {11:49:23 AM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
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
from fipy.meshes.numMesh.grid2D import Grid2D
import MA
from fipy.meshes.numMesh.testMeshBase import TestMeshBase
import fipy.tools.dump as dump

class TestGrid(TestMeshBase):
    def setUp(self):
	dx = 0.5
	dy = 2.
	nx = 3
	ny = 2
	
        self.mesh = Grid2D(nx = nx, ny = ny, dx = dx, dy = dy)     
        
        self.vertices = Numeric.array(((0., 0.), (1., 0.), (2., 0.), (3., 0.),
                                       (0., 1.), (1., 1.), (2., 1.), (3., 1.),
                                       (0., 2.), (1., 2.), (2., 2.), (3., 2.)))
        self.vertices = self.vertices * Numeric.array((dx, dy))
        self.faces = Numeric.array(((1, 0), (2, 1), (3, 2),
                                    (4, 5), (5, 6), (6, 7),
                                    (8, 9), (9, 10), (10, 11),
                                    (0, 4), (5, 1), (6, 2), (7, 3),
                                    (4, 8), (9, 5), (10, 6), (11, 7)))
        self.cells = Numeric.array(((0, 10, 3, 9),
				  (1 , 11, 4, 10),
				  (2, 12, 5, 11),
				  (3, 14, 6, 13),
				  (4, 15, 7, 14),
				  (5, 16, 8, 15)))


        self.externalFaces = Numeric.array((0, 1, 2, 6, 7, 8, 9 , 12, 13, 16))
        self.internalFaces = Numeric.array((3, 4, 5, 10, 11, 14, 15))
##        self.faceCellIds = Numeric.array(((0, 0), (1, 1), (2, 2),
##                                         (0, 3), (1, 4), (2, 5),
##                                         (3, 3), (4, 4), (5, 5),
##                                         (0, 0), (0, 1), (1, 2), (2, 2),
##                                         (3, 3), (3, 4), (4, 5), (5, 5)))
        self.faceCellIds = MA.masked_values(((0, -1), (1, -1), (2, -1),
                                          (0, 3), (1, 4), (2, 5),
                                          (3, -1), (4, -1), (5, -1),
                                          (0, -1), (0, 1), (1, 2), (2, -1),
                                          (3, -1), (3, 4), (4, 5), (5, -1)), -1)
        
        self.faceAreas = Numeric.array((dx, dx, dx, dx, dx, dx, dx, dx, dx,
                                        dy, dy, dy, dy, dy, dy, dy, dy))
        faceCoords = Numeric.take(self.vertices, self.faces)
        self.faceCenters = (faceCoords[:,0] + faceCoords[:,1]) / 2.

        self.faceNormals = Numeric.array(((0., -1.), (0., -1.), (0., -1.),
                                         (0., 1.), (0., 1.), (0., 1.),
                                         (0., 1.), (0., 1.), (0., 1.),
                                         (-1., 0), (1., 0), (1., 0), (1., 0),
                                         (-1., 0), (1., 0), (1., 0), (1., 0)))
        
	self.cellToFaceOrientations = Numeric.array(((1, 1, 1, 1), (1, 1, 1, -1), (1, 1, 1, -1),
						    (-1, 1, 1, 1), (-1, 1, 1, -1), (-1, 1, 1, -1)))
					 
	self.cellVolumes = Numeric.array((dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy))

	self.cellCenters = Numeric.array(((dx/2.,dy/2.), (3.*dx/2.,dy/2.), (5.*dx/2.,dy/2.),
                                          (dx/2.,3.*dy/2.), (3.*dx/2.,3.*dy/2.), (5.*dx/2.,3.*dy/2.)))
					  
	self.faceToCellDistances = MA.masked_values(((dy / 2., -1), (dy / 2., -1), (dy / 2., -1),
						    (dy / 2., dy / 2.), (dy / 2., dy / 2.), (dy / 2., dy / 2.),
						    (dy / 2., -1), (dy / 2., -1), (dy / 2., -1),
						    (dx / 2., -1), (dx / 2., dx / 2.), (dx / 2., dx / 2.),
						    (dx / 2., -1),
						    (dx / 2., -1), (dx / 2., dx / 2.), (dx / 2., dx / 2.),
						    (dx / 2., -1)), -1)
					  
	self.cellDistances = Numeric.array((dy / 2., dy / 2., dy / 2.,
					    dy, dy, dy,
					    dy / 2., dy / 2., dy / 2.,
					    dx / 2., dx, dx,
					    dx / 2.,
					    dx / 2., dx, dx,
					    dx / 2.))
        
        self.faceToCellDistanceRatios = self.faceToCellDistances[:,0] / self.cellDistances

        self.areaProjections = self.faceNormals * self.faceAreas[:,Numeric.NewAxis]

        self.facesBottom = (0, 1, 2)
	self.facesTop = (6, 7, 8)
	self.facesLeft = (9, 13)
	self.facesRight = (12, 16)

        self.tangents1 = Numeric.array(((1., 0), (1., 0),(1., 0),
                                       (-1., 0), (-1., 0),(-1., 0),
                                       (-1., 0), (-1., 0),(-1., 0),
                                       (0., -1.), (0., 1.), (0., 1.), (0., 1.),
                                       (0., -1.), (0., 1.), (0., 1.), (0., 1.)))

        self.tangents2 = Numeric.array(((0., 0), (0., 0),(0., 0),
                                       (-0., 0), (-0., 0),(-0., 0),
                                       (-0., 0), (-0., 0),(-0., 0),
                                       (0., -0.), (0., 0.), (0., 0.), (0., 0.),
                                       (0., -0.), (0., 0.), (0., 0.), (0., 0.)))

        self.cellToCellIDs = MA.masked_values(((-1, 1, 3, -1),
                                      (-1, 2, 4, 0),
                                      (-1, -1, 5, 1),
                                      (0, 4, -1, -1),
                                      (1, 5, -1, 3),
                                      (2, -1, -1, 4)), -1)

        self.cellToCellDistances = MA.masked_values(((dy / 2., dx, dy, dx / 2.),
                                                     (dy/ 2., dx, dy, dx),
                                                     (dy / 2., dx / 2., dy, dx),
                                                     (dy, dx, dy / 2., dx / 2.),
                                                     (dy, dx, dy / 2., dx),
                                                     (dy, dx / 2., dy / 2., dx)), -1)

        self.interiorCellIDs = Numeric.array(())

        self.exteriorCellIDs = Numeric.array((0, 1, 2, 3, 4, 5))

        self.cellNormals = Numeric.array( ( (  (0, -1), (1, 0), (0, 1), (-1, 0) ),
                                            (  (0, -1), (1, 0), (0, 1), (-1, 0) ),
                                            (  (0, -1), (1, 0), (0, 1), (-1, 0) ),
                                            (  (0, -1), (1, 0), (0, 1), (-1, 0) ),
                                            (  (0, -1), (1, 0), (0, 1), (-1, 0) ),
                                            (  (0, -1), (1, 0), (0, 1), (-1, 0) )  ) )

        vv = Numeric.array(((0, -dx), (dy, 0), (0, dx), (-dy, 0)))
        
        self.cellAreaProjections = Numeric.array(((vv,vv,vv,vv,vv,vv)))
                                       
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

class TestGridPickle(TestGrid):
    def setUp(self):
        TestGrid.setUp(self)
        import tempfile
        import os
	(f, fileName) = tempfile.mkstemp('.gz')
	pickledMesh = dump.write(self.mesh, fileName)
	self.mesh = dump.read(fileName)
        os.close(f)
	os.remove(fileName)

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestGrid))
    theSuite.addTest(unittest.makeSuite(TestGridPickle))
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')
