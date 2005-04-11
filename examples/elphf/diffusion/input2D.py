#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input2D.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 4/5/05 {8:10:16 PM} 
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

""" 
The same three-component diffusion problem as introduced in::
`examples/elphf/diffusion/input1D.py` but in 2D:
    
    >>> nx = 40
    >>> dx = 1.
    >>> L = nx * dx
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx = dx, dy = dx, nx = nx, ny = nx)
    
Again, the parameters are

    >>> parameters = {
    ...     'time step duration': 10000,
    ...     'solvent': {
    ...         'standard potential': 0.,
    ...         'barrier height': 0.
    ...     }
    ... }

    >>> parameters['substitutionals'] = (
    ...     {
    ...         'name': "c1",
    ...         'diffusivity': 1.,
    ...         'standard potential': 1.,
    ...         'barrier height': 1.
    ...     },
    ...     {
    ...         'name': "c2",
    ...         'diffusivity': 1.,
    ...         'standard potential': 1.,
    ...         'barrier height': 1.
    ...     }
    ... )
    
We again use ElPhF to create the variable fields

    >>> import fipy.models.elphf.elphf as elphf
    >>> fields = elphf.makeFields(mesh = mesh, 
    ...                           parameters = parameters)
    
and we separate the solution domain into two different concentration regimes
    
    >>> setCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > L/2)
    >>> fields['substitutionals'][0].setValue(0.3)
    >>> fields['substitutionals'][0].setValue(0.6,setCells)
    >>> fields['substitutionals'][1].setValue(0.6)
    >>> fields['substitutionals'][1].setValue(0.3,setCells)

We use ElPhF to create the governing equations for the fields

    >>> elphf.makeEquations(fields = fields, 
    ...                     parameters = parameters)
    
If we are running interactively, we create viewers to see the results 

    >>> if __name__ == '__main__':
    ...     import fipy.viewers
    ...     viewers = [fipy.viewers.make(vars = field) for field in fields['all']]
    ...     for viewer in viewers:
    ...         viewer.plot()

Now, we iterate the problem to equilibrium, plotting as we go

    >>> from fipy.solvers.linearLUSolver import LinearLUSolver
    >>> solver = LinearLUSolver()
    
    >>> for i in range(40):
    ...     for field in fields['substitutionals']:
    ...         field.updateOld()
    ...     for field in fields['substitutionals']:
    ...         field.equation.solve(var = field, 
    ...                              dt = parameters['time step duration'],
    ...                              solver = solver)
    ...     if __name__ == '__main__':
    ...         for viewer in viewers:
    ...             viewer.plot()

Since there is nothing to maintain the concentration separation in this problem, 
we verify that the concentrations have become uniform

    >>> print fields['substitutionals'][0].allclose(0.45, rtol = 1e-7, atol = 1e-7)
    1
    >>> print fields['substitutionals'][1].allclose(0.45, rtol = 1e-7, atol = 1e-7)
    1
"""
__docformat__ = 'restructuredtext'
 
## from fipy.tools.profiler.profiler import Profiler
## from fipy.tools.profiler.profiler import calibrate_profiler

## fudge = calibrate_profiler(10000)
## profile = Profiler('profile', fudge=fudge)

## profile.stop()

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    
    raw_input("finished")
    
