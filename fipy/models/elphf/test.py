#!/usr/bin/env python

## 
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/10/03 {3:23:47 PM}
 #                                last update: 12/29/03 {11:49:59 AM} 
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
 
from __future__ import nested_scopes

import Numeric

import unittest
from tests.testBase import TestBase
from meshes.grid2D import Grid2D
import elphf

class TestElPhF(TestBase):
    """
    Simple test case for the phase field equation.
    """
    def setUp(self):
	self.parameters = {
	    'diffusivity': 1.,
	    'time step duration': 10000.	
	}
	
	self.tolerance = 1e-7
	self.steps = 40
	
	self.getParameters(self.parameters)
	
	self.mesh = Grid2D(
	    dx = self.parameters['mesh']['dx'],
	    dy = self.parameters['mesh']['dy'],
	    nx = self.parameters['mesh']['nx'],
	    ny = self.parameters['mesh']['ny'])
	    
	fields = elphf.makeFields(self.mesh, self.parameters)
	    
	self.it = elphf.makeIterator(mesh = self.mesh, fields = fields, parameters = self.parameters)
	
    def testResult(self):
	self.it.iterate(steps = self.steps)

# 	for field in [self.parameters['phase']] + list(self.parameters['substitutionals']):
# 	    field['var'].plot()
# 	raw_input()
	    
	for field in [self.parameters['phase']] + list(self.parameters['substitutionals']):
	    final = field['var'].copy()
	    final.setValue(field['final'][0])
	    for fin in field['final'][1:]:
		setCells = self.mesh.getCells(fin['func'])
		final.setValue(fin['value'], setCells)

	    self.assertWithinTolerance(field['var'][0], final[0], self.tolerance)	
	    self.assertWithinTolerance(field['var'][-1], final[-1], self.tolerance)	
	    
class TestElPhF1D(TestElPhF):
    def getParameters(self, parameters):
	mesh = {
	    'nx': 40,
	    'ny': 1,
	    'dx': 1.,
	    'dy': 1.
	}
	
	L = mesh['nx'] * mesh['dx']
	
	parameters['mesh'] = mesh
	
	parameters['phase'] = {
	    'name': "xi",
	    'mobility': 1.,
	    'gradient energy': 1.,
	    'initial': (1.,),
	    'final': (1.,)
	}
	
	parameters['solvent'] = {
	    'standard potential': 0.,
	    'barrier height': 0.
	}
	
	rightFunc = lambda cell: cell.getCenter()[0] > L/2
	
	parameters['substitutionals'] = (
	    {
		'name': "c1",
		'standard potential': 1.,
		'barrier height': 1.,
		'initial': (
		    0.3,
		    {
			'value': 0.6,
			'func': rightFunc
		    }
		),
		'final': (0.45,)
	    },
	    {
		'name': "c2",
		'standard potential': 1.,
		'barrier height': 1.,
		'initial': (
		    0.6,
		    {
			'value': 0.3,
			'func': rightFunc
		    }
		),
		'final': (0.45,)
	    }
	)

class TestElPhF2D(TestElPhF):
    def getParameters(self, parameters):
	mesh = {
	    'nx': 40,
	    'ny': 40,
	    'dx': 1.,
	    'dy': 1.
	}
	
	L = mesh['nx'] * mesh['dx']
	
	parameters['mesh'] = mesh
	
	parameters['phase'] = {
	    'name': "xi",
	    'mobility': 1.,
	    'gradient energy': 1.,
	    'initial': (1.,),
	    'final': (1.,)
	}
	
	parameters['solvent'] = {
	    'standard potential': 0.,
	    'barrier height': 0.
	}
	
	rightFunc = lambda cell: cell.getCenter()[0] > L/2

	parameters['substitutionals'] = (
	    {
		'name': "c1",
		'standard potential': 1.,
		'barrier height': 1.,
		'initial': (
		    0.3,
		    {
			'value': 0.6,
			'func': rightFunc
		    }
		),
		'final': (0.45,)
	    },
	    {
		'name': "c2",
		'standard potential': 1.,
		'barrier height': 1.,
		'initial': (
		    0.6,
		    {
			'value': 0.3,
			'func': rightFunc
		    }
		),
		'final': (0.45,)
	    }
	)

class TestElPhF2DCorner(TestElPhF):
    def getParameters(self, parameters):
	mesh = {
	    'nx': 40,
	    'ny': 40,
	    'dx': 1.,
	    'dy': 1.
	}
	
	L = mesh['nx'] * mesh['dx']
	
	parameters['mesh'] = mesh
	
	parameters['phase'] = {
	    'name': "xi",
	    'mobility': 1.,
	    'gradient energy': 1.,
	    'initial': (1.,),
	    'final': (1.,)
	}
	
	parameters['solvent'] = {
	    'standard potential': 0.,
	    'barrier height': 0.
	}
	
	def cornerFunc(cell):
	    center = cell.getCenter()
	    return (center[0] > L/2) and (center[1] > L/2) 
	
	parameters['substitutionals'] = (
	    {
		'name': "c1",
		'standard potential': 1.,
		'barrier height': 1.,
		'initial': (
		    0.3,
		    {
			'value': 0.6,
			'func': cornerFunc
		    }
		),
		'final': (0.375,)
	    },
	    {
		'name': "c2",
		'standard potential': 1.,
		'barrier height': 1.,
		'initial': (
		    0.6,
		    {
			'value': 0.3,
			'func': cornerFunc
		    }
		),
		'final': (0.525,)
	    }
	)
	
class TestElPhF1Dphase(TestElPhF):
    def getParameters(self, parameters):
	self.tolerance = 1e-4
	self.steps = 10
	
	mesh = {
	    'nx': 400,
	    'ny': 1,
	    'dx': 0.01,
	    'dy': 0.01
	}
	
	L = mesh['nx'] * mesh['dx']
	
	parameters['mesh'] = mesh
	
	parameters['time step duration'] = 10000.	
	
	rightFunc = lambda cell: cell.getCenter()[0] > L/2
	
	parameters['phase'] = {
	    'name': "xi",
	    'mobility': 1.,
	    'gradient energy': 0.025,
	    'initial': (
		1.,
		{
		    'value': 0.,
		    'func': rightFunc
		}
	    )
	}
	
	parameters['solvent'] = {
	    'standard potential': 0.,
	    'barrier height': 1.
	}
	
	parameters['substitutionals'] = ()
	
    def testResult(self):
	self.it.iterate(steps = self.steps)

	field  = self.parameters['phase']
    
	x = Numeric.arange(float(self.parameters['mesh']['nx'])) 
	x -= (self.parameters['mesh']['nx'] - 1.) / 2.
	x *= self.parameters['mesh']['dx']
	d = Numeric.sqrt(field['gradient energy'] / (self.parameters['solvent']['barrier height']))
	final = (1. - Numeric.tanh(x/(2.*d))) / 2.
	
	self.assertArrayWithinTolerance(field['var'][:], final, self.tolerance)
	
class TestElPhF1DphaseBinary(TestElPhF):
    def getParameters(self, parameters):
	self.tolerance = 2e-3
	
	mesh = {
	    'nx': 400,
	    'ny': 1,
	    'dx': 0.01,
	    'dy': 0.01
	}
	
	L = mesh['nx'] * mesh['dx']
	
	parameters['mesh'] = mesh
	
	rightFunc = lambda cell: cell.getCenter()[0] > L/2
	
	parameters['phase'] = {
	    'name': "xi",
	    'mobility': 1.,
	    'gradient energy': 0.1,
	    'initial': (
		1.,
		{
		    'value': 0.,
		    'func': rightFunc
		}
	    ),
	    'final': (
		1.,
		{
		    'value': 0.,
		    'func': rightFunc
		}
	    ),
	}
	
	parameters['solvent'] = {
	    'standard potential': Numeric.log(.7/.3),
	    'barrier height': 1.
	}
	
	parameters['substitutionals'] = (
	    {
		'name': "c1",
		'standard potential': Numeric.log(.3/.7),
		'barrier height': parameters['solvent']['barrier height'], 
		'initial': (0.5,),
		'final': (
		    0.7,
		    {
			'value': 0.3,
			'func': rightFunc
		    }
		),
	    },
	)
	
class TestElPhF1DphaseQuaternary(TestElPhF):
    def getParameters(self, parameters):
	self.tolerance = 2e-3
	
	mesh = {
	    'nx': 400,
	    'ny': 1,
	    'dx': 0.01,
	    'dy': 0.01
	}
	
	L = mesh['nx'] * mesh['dx']
	
	parameters['mesh'] = mesh
	
	rightFunc = lambda cell: cell.getCenter()[0] > L/2
	
	parameters['phase'] = {
	    'name': "xi",
	    'mobility': 1.,
	    'gradient energy': 0.025,
	    'initial': (
		1.,
		{
		    'value': 0.,
		    'func': rightFunc
		}
	    ),
	    'final': (
		1.,
		{
		    'value': 0.,
		    'func': rightFunc
		}
	    ),
	}
	
	parameters['solvent'] = {
	    'standard potential': Numeric.log(.1/.2),
	    'barrier height': 1.
	}
	
	parameters['substitutionals'] = (
	    {
		'name': "c1",
		'standard potential': Numeric.log(.3/.4),
		'barrier height': parameters['solvent']['barrier height'],
		'initial': (0.35,),
		'final': (
		    0.4,
		    {
			'value': 0.3,
			'func': rightFunc
		    }
		)
	    },
	    {
		'name': "c2",
		'standard potential': Numeric.log(.4/.3),
		'barrier height': parameters['solvent']['barrier height'],
		'initial': (0.35,),
		'final': (
		    0.3,
		    {
			'value': 0.4,
			'func': rightFunc
		    }
		)
	    },
	    {
		'name': "c3",
		'standard potential': Numeric.log(.2/.1),
		'barrier height': parameters['solvent']['barrier height'],
		'initial': (0.15,),
		'final': (
		    0.1,
		    {
			'value': 0.2,
			'func': rightFunc
		    }
		)
	    }
	)
	    
class TestElPhF1DphaseTernaryAndElectrons(TestElPhF):
    def getParameters(self, parameters):
	self.tolerance = 2e-3
	
	mesh = {
	    'nx': 400,
	    'ny': 1,
	    'dx': 0.01,
	    'dy': 0.01
	}
	
	L = mesh['nx'] * mesh['dx']
	
	parameters['mesh'] = mesh
	
	rightFunc = lambda cell: cell.getCenter()[0] > L/2
	
	parameters['phase'] = {
	    'name': "xi",
	    'mobility': 1.,
	    'gradient energy': 0.025,
	    'initial': (
		1.,
		{
		    'value': 0.,
		    'func': rightFunc
		}
	    ),
	    'final': (
		1.,
		{
		    'value': 0.,
		    'func': rightFunc
		}
	    )
	}
	
	parameters['interstitials'] = (
	    {
		'name': "e-",
		'standard potential': Numeric.log(.3/.4),
		'barrier height': 0.,
		'initial': (0.35,),
		'final': (
		    0.4,
		    {
			'value': 0.3,
			'func': rightFunc
		    }
		)
	    },
	)
	    
	parameters['solvent'] = {
	    'standard potential': Numeric.log(.4/.6),
	    'barrier height': 1.
	}
	
	parameters['substitutionals'] = (
	    {
		'name': "c2",
		'standard potential': Numeric.log(.4/.3),
		'barrier height': parameters['solvent']['barrier height'],
		'initial': (0.35,),
		'final': (
		    0.3,
		    {
			'value': 0.4,
			'func': rightFunc
		    }
		)
	    },
	    {
		'name': "c3",
		'standard potential': Numeric.log(.2/.1),
		'barrier height': parameters['solvent']['barrier height'],
		'initial': (0.15,),
		'final': (
		    0.1,
		    {
			'value': 0.2,
			'func': rightFunc
		    }
		)
	    }
	)
	    
def suite():
    theSuite = unittest.TestSuite()
    theSuite.addTest(unittest.makeSuite(TestElPhF1D))
    theSuite.addTest(unittest.makeSuite(TestElPhF2D))
    theSuite.addTest(unittest.makeSuite(TestElPhF2DCorner))    
    theSuite.addTest(unittest.makeSuite(TestElPhF1Dphase))
    theSuite.addTest(unittest.makeSuite(TestElPhF1DphaseBinary))
    theSuite.addTest(unittest.makeSuite(TestElPhF1DphaseQuaternary))
    theSuite.addTest(unittest.makeSuite(TestElPhF1DphaseTernaryAndElectrons))
    return theSuite
    
if __name__ == '__main__':
    theSuite = suite()
    unittest.TextTestRunner(verbosity=2).run(theSuite)

            
            
