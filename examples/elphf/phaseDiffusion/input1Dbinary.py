#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 3/7/05 {5:16:26 PM} 
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
This example combines a phase field problem, as given in
``examples/elphf/input1Dphase.py``,
with a binary diffusion problem, such as described in the ternary example
``examples/elphf/input1D.py``,
on a 1D mesh

    >>> nx = 400
    >>> dx = 0.01
    >>> L = nx * dx
    >>> from fipy.meshes.grid1D import Grid1D
    >>> mesh = Grid1D(dx = dx, nx = nx)

The problem parameters are

    >>> parameters = {
    ...     'time step duration': 10000,
    ...     'substitutional molar volume': 1,
    ...     'phase': {
    ...         'name': "xi",
    ...         'mobility': 1.,
    ...         'gradient energy': 0.1,
    ...         'value': 1
    ...     },
    ... }

The thermodynamic parameters are chosen to give a solid phase rich 
in the solute and a liquid phase rich in the solvent

    >>> import Numeric
    >>> parameters['solvent'] = {
    ...     'standard potential': Numeric.log(.7/.3),
    ...     'barrier height': 1.
    ... }

    >>> parameters['substitutionals'] = (
    ...     {
    ...         'name': "c1",
    ...         'diffusivity': 1.,
    ...         'standard potential': Numeric.log(.3/.7),
    ...         'barrier height': parameters['solvent']['barrier height'], 
    ...     },
    ... )

We again let the ElPhF module create the appropriate fields and equations

    >>> import fipy.models.elphf.elphf as elphf
    >>> fields = elphf.makeFields(mesh = mesh, parameters = parameters)
    >>> elphf.makeEquations(fields = fields, parameters = parameters)

We start with a sharp phase boundary

.. raw:: latex

   $$ \xi =
   \begin{cases}
       1& \text{for $x \le L/2$,} \\
       0& \text{for $x > L/2$,}
   \end{cases} $$

or

    >>> setCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > L/2)
    >>> fields['phase'].setValue(1.)
    >>> fields['phase'].setValue(0.,setCells)

and with a uniform concentration field

.. raw:: latex

   $C_1 = 0.5$.
   
or

    >>> fields['substitutionals'][0].setValue(0.5)
    
If running interactively, we create viewers to display the results

    >>> if __name__ == '__main__':
    ...     import fipy.viewers
    ...
    ...     phaseViewer = fipy.viewers.make(vars = (fields['phase'],))
    ...     concViewer = fipy.viewers.make(vars = (fields['solvent'],) + fields['substitutionals'],
    ...                                    limits = {'datamin': 0, 'datamax': 1})
    ...     phaseViewer.plot()
    ...     concViewer.plot()

This problem does not have an analytical solution, so after iterating to
equilibrium

    >>> from fipy.solvers.linearLUSolver import LinearLUSolver
    >>> solver = LinearLUSolver()

    >>> for i in range(50):
    ...     for field in fields['all']:
    ...         field.updateOld()
    ...     fields['phase'].equation.solve(var = fields['phase'],
    ...                                    dt = parameters['time step duration'])
    ...     for field in fields['substitutionals']:
    ...         field.equation.solve(var = field, 
    ...                              dt = parameters['time step duration'],
    ...                              solver = solver)
    ...     if __name__ == '__main__':    
    ...         phaseViewer.plot()
    ...         concViewer.plot()

we confirm that the far-field phases have remained separated

    >>> ends = Numeric.take(fields['phase'], (0,-1))
    >>> Numeric.allclose(ends, (1.0, 0.0), rtol = 2e-3, atol = 2e-3)
    1
    
and that the concentration field has appropriately segregated into solute
rich and solute poor phases.

    >>> ends = Numeric.take(fields['substitutionals'][0], (0,-1))
    >>> Numeric.allclose(ends, (0.7, 0.3), rtol = 2e-3, atol = 2e-3)
    1
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    ## from fipy.tools.profiler.profiler import Profiler
    ## from fipy.tools.profiler.profiler import calibrate_profiler

    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)

    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
    
    # profile.stop()
	    
    raw_input("finished")

