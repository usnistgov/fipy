#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "test.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 4/1/05 {2:47:17 PM} 
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

from fipy.tests.doctestPlus import LateImportDocTestSuite
import fipy.tests.testProgram

def _suite():
    return LateImportDocTestSuite(docTestModuleNames = (
            'exponential1D.input',
            'exponential1D.tri2Dinput',
            'exponential2D.input',
            'exponential1DBack.input',
            'powerLaw1D.input',
            'exponential1DSource.input',
            'exponential2D.tri2Dinput',
            'exponential1DSource.tri2Dinput',
            'powerLaw1D.tri2Dinput'
        ), base = __name__)
    
if __name__ == '__main__':
    fipy.tests.testProgram.main(defaultTest='_suite')

            
            
