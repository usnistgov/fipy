#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:06:48 PM} 
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
from fipy.meshes.numMesh.mesh import Mesh
import MA
from mesh import Mesh
from fipy.meshes.numMesh.testMeshBase import TestMeshBase
import fipy.tools.dump as dump

class TestMesh3D(TestMeshBase):

    def setUp(self):
        dx = 2.
        dy = 1.23456
        dz = 1.e-1
        
        self.vertices = Numeric.array(((0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
                                       (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
                                       (2., 0., 0.), (2., 0., 1.)))
        
        self.faces = MA.masked_values(((0, 1, 2, 3), (7, 6, 5, 4),
                                    (3, 7, 4, 0), (5, 6, 2, 1),
                                    (1, 0, 4, 5), (3, 2, 6, 7),
                                    (1, 8, 2, -1), (9, 5, 6, -1), (8, 1, 5, 9), (8, 9, 6, 2)),-1)
        
        self.cells = MA.masked_values(((0, 1, 2, 3, 4, 5),
                                       (3 , 6, 7, 8, 9, -1)), -1)
    

        self.vertices = self.vertices * Numeric.array((dx, dy, dz))

        self.mesh = Mesh(self.vertices, self.faces, self.cells)
        
        self.externalFaces = Numeric.array((0, 1, 2, 4, 5, 6, 7, 8, 9))
        self.internalFaces = Numeric.array((3,))
        self.faceCellIds = MA.masked_values(((0, -1), (0, -1), (0, -1),
                                         (0, 1), (0, -1), (0, -1),
                                         (1, -1), (1, -1), (1, -1), (1, -1)), -1)

        dxdy = dx * dy
        dxdz = dx * dz
        dydz = dy * dz
        self.faceAreas = Numeric.array((dxdy, dxdy, dydz, dydz, dxdz, dxdz,
                                        dxdy/2., dxdy/2., dxdz, Numeric.sqrt(dx**2 + dy**2) * dz))
        faceCoords = Numeric.take(self.vertices, MA.filled(self.faces, 0))
        self.faceCenters = faceCoords[:,0] + faceCoords[:,1] + faceCoords[:,2] + faceCoords[:,3]
        numVex = Numeric.array((4., 4., 4., 4., 4., 4., 3., 3., 4., 4.))
        self.faceCenters[:,0] = self.faceCenters[:,0] / numVex
        self.faceCenters[:,1] = self.faceCenters[:,1] / numVex
        self.faceCenters[:,2] = self.faceCenters[:,2] / numVex

        self.faceNormals = Numeric.array(((0., 0., -1.),
                                          (0., 0., 1.),
                                          (-1, 0., 0.),
                                          (1. , 0., 0.),
                                          (0, -1., 0.),
                                          (0, 1., 0.),
                                          (0., 0., -1.),
                                          (0., 0., 1.),
                                          (0., -1., 0.),
                                          (dy / Numeric.sqrt(dy**2 + dx**2), dx / Numeric.sqrt(dy**2 + dx**2), 0.)))

        self.cellToFaceOrientations = MA.masked_values(((1, 1, 1, 1, 1, 1),
                                                         (-1, 1, 1, 1, 1, -2)), -2)
							 
	self.cellVolumes = Numeric.array((dx*dy*dz, dx*dy*dz / 2.))

	self.cellCenters = Numeric.array(((dx/2.,dy/2.,dz/2.), (dx+dx/3.,dy/3.,dz/2.)))

        d1 = Numeric.sqrt((dx / 3.)**2 + (dy / 6.)**2)
        d2 = Numeric.sqrt((dx / 6.)**2 + (dy / 3.)**2)
        d3 = Numeric.sqrt((dx / 6.)**2 + (dy / 6.)**2)
        d4 = Numeric.sqrt((5 * dx / 6.)**2 + (dy / 6.)**2)

        self.faceToCellDistances = MA.masked_values(((dz / 2., -1),
                                                    (dz / 2., -1),
                                                    (dx / 2., -1),
                                                    (dx / 2., d1),
                                                    (dy / 2., -1),
                                                    (dy / 2., -1),
                                                    (dz / 2., -1),
                                                    (dz / 2., -1),
                                                    (d2, -1),
                                                    (d3, -1)), -1)


                                                    
	self.cellDistances = Numeric.array((dz / 2., dz / 2., dx / 2.,
                                            d4,
					    dy / 2., dy / 2., dz / 2., dz / 2.,
					    d2,
					    d3))

        self.faceToCellDistanceRatios = self.faceToCellDistances[:,0] / self.cellDistances

        self.areaProjections = self.faceNormals * self.faceAreas[:,Numeric.NewAxis]
    
        v1 = Numeric.take(self.vertices, self.faces[:,0])
        tmp = self.faceCenters - v1
        self.tangents1 = tmp / fipy.tools.array.sqrtDot(tmp, tmp)[:,Numeric.NewAxis]
        
        tmp = fipy.tools.array.crossProd(self.tangents1, self.faceNormals)
        self.tangents2 = tmp / fipy.tools.array.sqrtDot(tmp, tmp)[:,Numeric.NewAxis]

        self.cellToCellIDs = MA.masked_values(((-1, -1, -1, 1, -1, -1),
                                               (0, -1, -1, -1, -1, -1)), -1)

        self.cellToCellDistances = MA.masked_values(((dz / 2., dz / 2., dx / 2., d4, dy / 2., dy / 2.),
                                                      (d4, -1, -1, -1, -1, -1)), -1)

class TestMesh3DPickle(TestMesh3D):
    def setUp(self):
        TestMesh3D.setUp(self)
        pickledMesh = dump.write(self.mesh, 'pickledMesh')
        self.mesh = dump.read('pickledMesh')

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestMesh3D))
    theSuite.addTest(unittest.makeSuite(TestMesh3DPickle))
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')
