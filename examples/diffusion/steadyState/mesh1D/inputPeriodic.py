#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 4/6/05 {4:32:13 PM} 
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

r"""

    >>> import examples.diffusion.steadyState.mesh1D.input
    >>> exec(examples.diffusion.steadyState.mesh1D.input.script())
    
One can then solve the same problem as in
`examples/diffusion/steadyState/mesh1D/input.py` but with periodic
boundary conditions.

   >>> from fipy.boundaryConditions.periodicBoundaryCondition import PeriodicBoundaryCondition
   >>> boundaryConditions = (PeriodicBoundaryCondition(mesh.getFacesLeft(), mesh.getFacesRight()),)

A `TransientTerm` is used to provide some fixed point, otherwise the
solver has no fixed value and can become unstable.
    
   >>> from fipy.terms.transientTerm import TransientTerm
   >>> eq = TransientTerm(coeff = 1e-7) - ImplicitDiffusionTerm()
   >>> eq.solve(var = var, boundaryConditions = boundaryConditions)

   >>> if __name__ == '__main__':
   ...     viewer.plot()
   ...     raw_input("press key to continue")

The result of the calculation will be the average value over the domain.

   >>> var.allclose((valueLeft + valueRight) / 2.).getValue()
   1
   
"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
    raw_input("finished")
