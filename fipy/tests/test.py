#!/usr/bin/env python

## 
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 11/26/03 {3:23:47 PM}
 #                                last update: 1/16/04 {11:57:01 AM} { 2:26:30 PM}
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

"""Run all the test cases
"""

import unittest

import fivol.examples.diffusion.variable.test
import fivol.examples.diffusion.steadyState.test
import fivol.examples.diffusion.explicit.test
import fivol.examples.convection.test
import fivol.examples.elphf.test
import fivol.examples.phase.examples.test

if __name__ == '__main__':
    theSuite = unittest.TestSuite()
    
    theSuite.addTest(fivol.examples.diffusion.steadyState.test.suite())
    theSuite.addTest(fivol.examples.diffusion.explicit.test.suite())
    theSuite.addTest(fivol.examples.diffusion.variable.test.suite())
    theSuite.addTest(fivol.examples.convection.test.suite())
    theSuite.addTest(fivol.examples.phase.examples.test.suite())
    theSuite.addTest(fivol.examples.elphf.test.suite())
    
    unittest.TextTestRunner(verbosity=2).run(theSuite)

