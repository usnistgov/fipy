#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input1Ddimensional.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 4/5/05 {8:09:56 PM} 
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
In this example, we present the same three-component diffusion problem 
introduced in ``examples/elphf/input1D.py``
but we demonstrate FiPy's facility to use dimensional quantities.

    >>> from fipy.tools.dimensions.physicalField import PhysicalField

We solve the problem on a 40 mm long 1D mesh

    >>> nx = 40
    >>> dx = PhysicalField(1.,"mm")
    >>> L = nx * dx
    >>> from fipy.meshes.grid1D import Grid1D
    >>> mesh = Grid1D(dx = dx, nx = nx)

The dimensional parameters for this problem are

    >>> parameters = {
    ...     'time step duration': "1000 s",
    ...     'solvent': {
    ...         'standard potential': 0.,
    ...         'barrier height': 0.
    ...     }
    ... }

    >>> parameters['substitutionals'] = (
    ...     {
    ...         'name': "c1",
    ...         'diffusivity': "1.e-9 m**2/s",
    ...         'standard potential': 1.,
    ...         'barrier height': 1.
    ...     },
    ...     {
    ...         'name': "c2",
    ...         'diffusivity': "1.e-9 m**2/s",
    ...         'standard potential': 1.,
    ...         'barrier height': 1.
    ...     }
    ... )

We use ElPhF to create the variable fields

    >>> import fipy.models.elphf.elphf as elphf
    >>> fields = elphf.makeFields(mesh = mesh, parameters = parameters)
    
and we separate the solution domain into two different concentration regimes
    
    >>> setCells = mesh.getCells(filter = lambda cell: 
    ...                          cell.getCenter()[0] > mesh.getPhysicalShape()[0]/2)
    >>> fields['substitutionals'][0].setValue("0.3 mol/m**3")
    >>> fields['substitutionals'][0].setValue("0.6 mol/m**3",setCells)
    >>> fields['substitutionals'][1].setValue("0.6 mol/m**3")
    >>> fields['substitutionals'][1].setValue("0.3 mol/m**3",setCells)

We use ElPhF again to create the governing equations for the fields

    >>> elphf.makeEquations(fields = fields, 
    ...                     parameters = parameters)
    
If we are running interactively, we create a viewer to see the results 

    >>> if __name__ == '__main__':
    ...     import fipy.viewers
    ...     viewer = fipy.viewers.make(
    ...         vars = (fields['solvent'],) + fields['substitutionals'],
    ...         limits = {'datamin': 0, 'datamax': 1})
    ...     viewer.plot()

Now, we iterate the problem to equilibrium, plotting as we go

    >>> from fipy.solvers.linearLUSolver import LinearLUSolver
    >>> solver = LinearLUSolver()

    >>> for i in range(40):
    ...     for field in fields['substitutionals']:
    ...         field.updateOld()
    ...     for field in fields['substitutionals']:
    ...         field.equation.solve(var = field, 
    ...                              dt = 10000, # (scaling doesn't work) parameters['time step duration'],
    ...                              solver = solver)
    ...     if __name__ == '__main__':
    ...         viewer.plot()

Since there is nothing to maintain the concentration separation in this problem, 
we verify that the concentrations have become uniform

    >>> print fields['substitutionals'][0].getScaled().allclose("0.45 mol/m**3",
    ...     atol = "1e-7 mol/m**3", rtol = 1e-7)
    1
    >>> print fields['substitutionals'][1].getScaled().allclose("0.45 mol/m**3",
    ...     atol = "1e-7 mol/m**3", rtol = 1e-7)
    1
    
.. note::
    
   The absolute tolerance `atol` must be in units compatible with the value to 
   be checked, but the relative tolerance `rtol` is dimensionless.
"""
__docformat__ = 'restructuredtext'


if __name__ == '__main__':
    ## from fipy.tools.profiler.profiler import Profiler
    ## from fipy.tools.profiler.profiler import calibrate_profiler

    # fudge = calibrate_profiler(10000)
    # profile = Profiler('profile', fudge=fudge)

    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    
    # profile.stop()
	    
    raw_input("finished")

