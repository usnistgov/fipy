#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 4/2/04 {4:05:48 PM} 
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
import doctest
import fipy.tests.testProgram
import fipy.meshes.numMesh.testGrid
import fipy.meshes.numMesh.testMesh
import fipy.meshes.numMesh.testMesh3D
import fipy.meshes.numMesh.testTri2D
import fipy.meshes.numMesh.mesh
import fipy.meshes.numMesh.gmshImport
import fipy.meshes.numMesh.adaptiveMesh
def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(fipy.meshes.numMesh.testGrid.suite())
    theSuite.addTest(fipy.meshes.numMesh.testMesh.suite())
    theSuite.addTest(fipy.meshes.numMesh.testMesh3D.suite())
    theSuite.addTest(fipy.meshes.numMesh.testTri2D.suite())
    theSuite.addTest(doctest.DocTestSuite(fipy.meshes.numMesh.mesh))
    theSuite.addTest(doctest.DocTestSuite(fipy.meshes.numMesh.gmshImport))
    theSuite.addTest(doctest.DocTestSuite(fipy.meshes.numMesh.adaptiveMesh))
    return theSuite
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='suite')
