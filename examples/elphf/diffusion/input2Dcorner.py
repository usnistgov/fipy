#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input2Dcorner.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 12/10/04 {5:15:24 PM} 
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
The same three-component diffusion problem introduced in
`examples/elphf/diffusion/input2D.py`, but with the initial step in
concentration occupying only one corner instead of half the domain.

We create a mesh

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

    >>> def cornerFunc(cell):
    ...     center = cell.getCenter()
    ...     return (center[0] > L/2) and (center[1] > L/2) 
        
    >>> setCells = mesh.getCells(filter = cornerFunc)
    >>> fields['substitutionals'][0].setValue(0.3)
    >>> fields['substitutionals'][0].setValue(0.6,setCells)
    >>> fields['substitutionals'][1].setValue(0.6)
    >>> fields['substitutionals'][1].setValue(0.3,setCells)

We use ElPhF to create the governing equations for the fields

    >>> elphf.makeEquations(fields = fields, 
    ...                     parameters = parameters)

If we are running interactively, we create viewers to see the results 

    >>> if __name__ == '__main__':
    ...     from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
    ...     viewers = [Grid2DGistViewer(var = field) for field in fields['all']]
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

    >>> fields['substitutionals'][0].allclose(0.375, rtol = 1e-7, atol = 1e-7)
    1
    >>> fields['substitutionals'][1].allclose(0.525, rtol = 1e-7, atol = 1e-7)
    1
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
    
    raw_input("finished")