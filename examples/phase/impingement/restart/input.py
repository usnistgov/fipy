#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 10/27/04 {9:52:33 AM}
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
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

r"""
In this example we solve the same set of equations as in::
    
    $ examples/phase/impingement/mesh20x20/input.py
    
but a restart method is demonstrated. 

We again initialize the system by running the base script
    
    >>> import examples.phase.impingement.mesh20x20.base
    >>> exec(examples.phase.impingement.mesh20x20.base.script())

but this time we only iterate for half as many time steps

    >>> for i in range(steps/2):
    ...     it.timestep(dt = timeStepDuration)
    ...     if __name__ == '__main__':
    ...         phaseViewer.plot()
    ...         thetaProductViewer.plot()

We confirm that the solution has not yet converged to that given by 
Ryo Kobayshi's FORTRAN code:

    >>> theta.allclose(testData, rtol = 1e-10, atol = 1e-10)
    0

We save the variables to disk

    >>> import fipy.tools.dump as dump
    >>> import tempfile
    >>> import os
    >>> tmp = tempfile.gettempdir()
    >>> fileName = os.path.join(tmp, 'data')
    >>> dump.write({'phase' : phase, 'theta' : theta, 'mesh' : mesh}, fileName)
   
and then recall them to test the data pickling mechanism

    >>> data = dump.read(fileName)
    >>> newPhase = data['phase']
    >>> newTheta = data['theta']
    >>> newMesh = data['mesh']
   
We rebuild the equations:

    >>> newThetaEq = ThetaEquation(
    ...     var = newTheta,
    ...     solver = LinearPCGSolver(
    ...         tolerance = 1.e-15, 
    ...         steps = 2000
    ...     ),
    ...     boundaryConditions = (
    ...         FixedFlux(newMesh.getExteriorFaces(), 0.),
    ...     ),
    ...     parameters = thetaParameters,
    ...     fields = {
    ...         'phase' : newPhase
    ...     }
    ... )
     
    >>> newPhaseEq = PhaseEquation(
    ...     var = newPhase,
    ...     mPhi = Type1MPhiVariable,
    ...     solver = LinearPCGSolver(
    ...         tolerance = 1.e-15, 
    ...         steps = 1000
    ...     ),
    ...     boundaryConditions = (
    ...         FixedFlux(newMesh.getExteriorFaces(), 0.),
    ...      ),
    ...      parameters = phaseParameters,
    ...      fields = {
    ...          'theta' : newTheta,
    ...          'temperature' : temperature
    ...      }
    ... )

the iterator:

    >>> newIt = Iterator((newThetaEq, newPhaseEq))
   
and, if the example is run interactively, we recreate the viewers:

    >>> if __name__ == '__main__':
    ...     from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
    ...     newPhaseViewer = Grid2DGistViewer(var = newPhase, palette = 'rainbow.gp', 
    ...                                       minVal = 0., maxVal = 1., grid = 0)
    ...     from Numeric import pi
    ...     newThetaProd = -pi + phase * (newTheta + pi)
    ...     newThetaProductViewer = \
    ...         Grid2DGistViewer(var = newThetaProd , palette = 'rainbow.gp', 
    ...                          minVal = -pi, maxVal = pi, grid = 0)
    ...     newPhaseViewer.plot()
    ...     newThetaProductViewer.plot()

and finish doing the iterations

    >>> for i in range(steps - steps/2):
    ...     newIt.timestep(dt = timeStepDuration)
    ...     if __name__ == '__main__':
    ...         newPhaseViewer.plot()
    ... 	newThetaProductViewer.plot()

Finally, we check that the results have at last converged to Ryo Kobayshi's 
FORTRAN code:

    >>> newTheta.allclose(testData, rtol = 1e-10, atol = 1e-10)
    1
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
    
    raw_input("finished")
