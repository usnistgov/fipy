#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 10/26/04 {8:45:59 AM} 
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
A simple 1D phase-field problem to test the `PhaseEquation` element of
ElPhF.

The single-component phase field governing equation can be represented as

.. raw:: latex

   $$ \frac{1}{M_\xi} \frac{\partial \xi}{\partial t} 
   =  \kappa_\xi \nabla^2 \xi - 2\xi(1-\xi)(1-2\xi) W $$

where 

.. raw:: latex

   $\xi$ is the phase field,
   $t$  is time,
   $M_\xi$ is the phase field mobility,
   $\kappa_\xi$ is the phase field gradient energy coefficient, and
   $W$ is the phase field barrier energy.
   
We solve the problem on a 1D mesh

    >>> nx = 400
    >>> dx = 0.01
    >>> ny = 1
    >>> dy = dx
    >>> L = nx * dx
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx = dx, dy = dy, nx = nx, ny = ny)

Rather than rewriting the same code in every electrochemistry example, 
we use the ElPhF module

    >>> import fipy.models.elphf.elphf as elphf

to build the approriate variable fields from

    >>> parameters = {
    ...     'time step duration': 10000,
    ...     'phase': {
    ... 	    'name': "xi",
    ... 	    'mobility': 1.,
    ... 	    'gradient energy': 0.025,
    ... 	    'value': 1.
    ...     },
    ...     'solvent': {
    ... 	    'standard potential': 0.,
    ... 	    'barrier height': 1.
    ...     }
    ... }
    
    >>> fields = elphf.makeFields(mesh = mesh, parameters = parameters)

We separate the phase field into electrode and electrolyte regimes

    >>> setCells = mesh.getCells(filter = lambda cell: cell.getCenter()[0] > L/2)
    >>> fields['phase'].setValue(1.)
    >>> fields['phase'].setValue(0.,setCells)

We use the ElPhF module again to create governing equations from the fields

    >>> equations = elphf.makeEquations(mesh = mesh, 
    ...                                 fields = fields, 
    ...                                 parameters = parameters)

If we are running interactively, we will want to see the results

    >>> if __name__ == '__main__':
    ...     from fipy.viewers.gist1DViewer import Gist1DViewer
    ...     viewer = Gist1DViewer(vars = (fields['phase'],))
    ...     viewer.plot()

Now, we iterate to equilibrium, plotting as we go

    >>> from fipy.iterators.iterator import Iterator
    >>> it = Iterator(equations = equations)
    >>> for i in range(50):
    ...     it.timestep(1)
    ...     if __name__ == '__main__':
    ...         viewer.plot()

The phase field has the expected analytical form

.. raw:: latex

   $$ \xi(x) = \frac{1}{2}(1 - \tanh\frac{x - L/2}{2d}) $$
   
where the interfacial thickness is given by

.. raw:: latex

   $ d = \sqrt{\kappa_{\xi}/W} $.
   
We verify that the correct equilibrium solution is attained

    >>> x = mesh.getCellCenters()[:,0]
    
    >>> import Numeric
    >>> d = Numeric.sqrt(parameters['phase']['gradient energy']
    ...     / (parameters['solvent']['barrier height']))
    >>> analyticalArray = (1. - Numeric.tanh((x - L/2.)/(2.*d))) / 2.

    >>> fields['phase'].allclose(analyticalArray, rtol = 1e-4, atol = 1e-4)
    1
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
##     from fipy.tools.profiler.profiler import Profiler
##     from fipy.tools.profiler.profiler import calibrate_profiler

    import doctest
    doctest.testmod()

##     fudge = calibrate_profiler(10000)
##     profile = Profiler('profile', fudge=fudge)

##     profile.stop()
	    
    raw_input("finished")

