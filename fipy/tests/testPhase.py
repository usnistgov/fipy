#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "testSteadyStateDiffusion.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 11/28/03 {10:46:22 AM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
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
 # ========================================================================
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

"""Test steady-state diffusion solutions
"""
 
import unittest
import os
import cPickle
import tests
from meshes.grid2D import Grid2D
from phase.phaseEquation import PhaseEquation
from solvers.linearPCGSolver import LinearPCGSolver
from boundaryConditions.fixedValue import FixedValue
from boundaryConditions.fixedFlux import FixedFlux
from iterators.iterator import Iterator
from variables.variable import Variable

class TestPhase(unittest.TestCase):
    """
    Simple test case for the phase field equation.
    """
    def setUp(self):

        phaseParameters={
            'tau' :        0.1,
            'epsilon' :    0.008,
            's' :          0.01,
            'alpha' :      0.015,
            'c2':          0.0,
            'anisotropy':  0.,
            'symmetry':    4.
            }
        
        valueLeft=1.
        valueRight=1.
        
        L = 1.5
        nx = 100
        dx = L/nx
        
        mesh = Grid2D(dx,1.,nx,1)
        
        self.phase = Variable(
            name = 'PhaseField',
            mesh = mesh,
            value = 1.
            )
        
        theta = Variable(
            name = 'Theta',
            mesh = mesh,
            value = 1.,
            hasOld = 0
            )
        
        def func(x,L=L):
            if x[0] > L / 2.:
                return 1
            else:
                return 0

        rightCells = mesh.getCells(func)
        
        theta.setValue(0.,rightCells)

        eq = PhaseEquation(
            self.phase,
            theta = theta,
            temperature = 1.,
            solver = LinearPCGSolver(
            tolerance = 1.e-15, 
            steps = 1000
            ),
            boundaryConditions=(
            FixedValue(mesh.getFacesLeft(),valueLeft),
            FixedValue(mesh.getFacesRight(),valueRight)),
            parameters = phaseParameters
            )
        
        self.it = Iterator((eq,))

    def assertWithinTolerance(self, first, second, tol = 1e-10, msg=None):
        """Fail if the two objects are unequal by more than tol.
        """
        if abs(first - second) > tol:
            raise self.failureException, (msg or '%s !~ %s' % (first, second))
        
    def testResult(self):
        self.it.iterate(100,0.02)
        array = self.phase.getArray()

        filestream=os.popen('gunzip --fast -c < %s/testPhaseData.gz'%tests.__path__[0],'r')
        testArray=cPickle.load(filestream)
        filestream.close()

        for i in range(len(array)):
            norm = abs(array[i] - testArray[0,i])        
            self.assertWithinTolerance(norm, 0.0, 1e-8,("cell(%g)'s value of %g differs from %g by %g" % (i,array[i],testArray[0,i],norm)))
	
def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestPhase))
    return theSuite
    
if __name__ == '__main__':
    theSuite = suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)

            
            
