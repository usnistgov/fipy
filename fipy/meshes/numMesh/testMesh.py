#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 3/5/04 {11:14:15 AM} 
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
import fivol.tests.testProgram
import Numeric
import MA
from fivol.meshes.mesh2D import Mesh2D
from fivol.meshes.testMeshBase import TestMeshBase
import fivol.tools.dump as dump

class TestMesh(TestMeshBase):
    def setUp(self):
        dx = 0.5
        dy = 2.
        nx = 3
        ny = 2
        self.vertices = Numeric.array(((0., 0.), (1., 0.), (2., 0.), (3., 0.),
                                       (0., 1.), (1., 1.), (2., 1.), (3., 1.),
                                       (0., 2.), (1., 2.), (2., 2.), (3., 2.),
                                       (4., 1.)))
        self.vertices = self.vertices * Numeric.array((dx, dy))
        self.faces = Numeric.array(((1, 0), (2, 1), (3, 2),
                                    (4, 5), (5, 6), (6, 7),
                                    (8, 9), (9, 10), (10, 11),
                                    (0, 4), (5, 1), (6, 2), (7, 3),
                                    (4, 8), (9, 5), (10, 6), (11, 7),
                                    (12, 3), (7, 12), (11, 12)))
        
        self.cells = MA.masked_values(((0, 10, 3, 9),
				  (1 , 11, 4, 10),
				  (2, 12, 5, 11),
				  (3, 14, 6, 13),
				  (4, 15, 7, 14),
				  (5, 16, 8, 15),
                                  (17, 18, 12, -1),
                                  (18, 19 ,16, -1)), -1)

        self.mesh = Mesh2D(self.vertices, self.faces, self.cells)

        self.externalFaces = Numeric.array((0, 1, 2, 6, 7, 8, 9, 13, 17, 19))
        self.internalFaces = Numeric.array((3, 4, 5, 10, 11, 12, 14, 15, 16, 18))
        self.faceCellIds = MA.masked_values(((0, -1), (1, -1), (2, -1),
                                         (0, 3), (1, 4), (2, 5),
                                         (3, -1), (4, -1), (5, -1),
                                         (0, -1), (0, 1), (1, 2), (2, 6),
                                         (3, -1), (3, 4), (4, 5), (5, 7),
                                         (6, -1), (6, 7), (7, -1)), -1)
        self.faceAreas = Numeric.array((dx, dx, dx, dx, dx, dx, dx, dx, dx,
                                       dy, dy, dy, dy, dy, dy, dy, dy,
                                       Numeric.sqrt(dx**2 + dy**2), dx, Numeric.sqrt(dx**2 + dy**2)))
        faceCoords = Numeric.take(self.vertices, self.faces)
        self.faceCenters = (faceCoords[:,0] + faceCoords[:,1]) / 2.

        self.faceNormals = Numeric.array(((0., -1.), (0., -1.), (0., -1.),
                                         (0., 1.), (0., 1.), (0., 1.),
                                         (0., 1.), (0., 1.), (0., 1.),
                                         (-1., 0), (1., 0), (1., 0), (1., 0),
                                         (-1., 0), (1., 0), (1., 0), (1., 0),
                                         (dy / Numeric.sqrt(dx**2 + dy**2), -dx / Numeric.sqrt(dx**2 + dy**2)),
                                         (0., 1.),
                                         (dy / Numeric.sqrt(dx**2 + dy**2), dx / Numeric.sqrt(dx**2 + dy**2))))

	self.cellToFaceOrientations = MA.masked_values(((1, 1, 1, 1), (1, 1, 1, -1), (1, 1, 1, -1),
							(-1, 1, 1, 1), (-1, 1, 1, -1), (-1, 1, 1, -1),
							(1, 1, -1, 0), (-1, 1, -1, 0)), 0)

	self.cellVolumes = Numeric.array((dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy / 2., dx*dy / 2.))

        self.cellCenters = Numeric.array(((dx/2.,dy/2.), (3.*dx/2.,dy/2.), (5.*dx/2.,dy/2.), 
                                          (dx/2.,3.*dy/2.), (3.*dx/2.,3.*dy/2.), (5.*dx/2.,3.*dy/2.),
                                          (3.*dx+dx/3.,2.*dy/3.), (3.*dx+dx/3.,4.*dy/3.)))

	self.faceToCellDistances = MA.masked_values(((dy / 2., -1), (dy / 2., -1), (dy / 2., -1),
						    (dy / 2., dy / 2.), (dy / 2., dy / 2.), (dy / 2., dy / 2.),
						    (dy / 2., -1), (dy / 2., -1), (dy / 2., -1),
						    (dx / 2., -1), (dx / 2., dx / 2.), (dx / 2., dx / 2.),
						    (dx / 2., Numeric.sqrt((dx / 3.)**2 + (dy / 6.)**2)),
						    (dx / 2., -1), (dx / 2., dx / 2.), (dx / 2., dx / 2.),
						    (dx / 2., Numeric.sqrt((dx / 3.)**2 + (dy / 6.)**2)),
						    (Numeric.sqrt((dx / 6.)**2 + (dy / 6.)**2), -1),
						    (Numeric.sqrt((dx / 6.)**2 + (dy / 3.)**2), Numeric.sqrt((dx / 6.)**2 + (dy / 3.)**2)),
						    (Numeric.sqrt((dx / 6.)**2 + (dy / 6.)**2), -1)), -1)
						    
	self.cellDistances = Numeric.array((dy / 2., dy / 2., dy / 2.,
						dy, dy, dy,
						dy / 2., dy / 2., dy / 2.,
						dx / 2., dx, dx,
						Numeric.sqrt((dx / 3. + dx / 2.)**2 + (dy / 6.)**2),
						dx / 2., dx, dx,
						Numeric.sqrt((dx / 3. + dx / 2.)**2 + (dy / 6.)**2),
						Numeric.sqrt((dx / 6.)**2 + (dy / 6.)**2),
						2. * dy / 3.,
						Numeric.sqrt((dx / 6.)**2 + (dy / 6.)**2)))
						    

        self.faceToCellDistanceRatios = self.faceToCellDistances[:,0] / self.cellDistances

        self.areaProjections = self.faceNormals * self.faceAreas[:,Numeric.NewAxis]

        self.tangents1 = Numeric.array(((1., 0), (1., 0),(1., 0),
                                        (-1., 0), (-1., 0),(-1., 0),
                                        (-1., 0), (-1., 0),(-1., 0),
                                        (0., -1.), (0., 1.), (0., 1.), (0., 1.),
                                        (0., -1.), (0., 1.), (0., 1.), (0., 1.),
                                        (dx / Numeric.sqrt(dx**2 +dy**2), dy / Numeric.sqrt(dx**2 +dy**2)),
                                        (-1, 0.),
                                        (-dx / Numeric.sqrt(dx**2 +dy**2), dy / Numeric.sqrt(dx**2 +dy**2))))
        

        self.tangents2 = Numeric.array(((0., 0), (0., 0),(0., 0),
                                        (-0., 0), (-0., 0),(-0., 0),
                                        (-0., 0), (-0., 0),(-0., 0),
                                        (0., -0.), (0., 0.), (0., 0.), (0., 0.),
                                        (0., -0.), (0., 0.), (0., 0.), (0., 0.),
                                        (0., 0), (0., 0),(0., 0)))


class TestMeshPickle(TestMesh):
    def setUp(self):
        TestMesh.setUp(self)
        pickledMesh = dump.write(self.mesh, 'pickledMesh')
##        self.mesh = dump.read('pickledMesh')

def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestMesh))
    theSuite.addTest(unittest.makeSuite(TestMeshPickle))
    return theSuite
    
if __name__ == '__main__':
    fivol.tests.testProgram.main(defaultTest='suite')
