#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 3/7/05 {2:53:22 PM} { 1:23:41 PM}
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

Here we create a level set variable in one dimension. The level set
variable calculates its value so that that value at any point in the
domain is the distance from the zero level set. This can be
represented succinctly in the following equation with a boundary
condition at the zero level set such that,

.. raw:: latex

    $$ \frac{\partial \phi}{\partial x} = 1 $$

with the boundary condition,

.. raw:: latex

    $\phi = 0$ at $x = L / 2$.

The solution to this problem will be demonstrated in the following
script. Firstly, setup the parameters.

   >>> dx = 0.5
   >>> dy = 2.
   >>> nx = 10
   >>> ny = 1

Construct the mesh.

   >>> from fipy.meshes.grid1D import Grid1D
   >>> mesh = Grid1D(dx = dx, nx = nx)

Construct a `distanceVariable` object.

   >>> from fipy.models.levelSet.distanceFunction.distanceVariable import DistanceVariable
   >>> var = DistanceVariable(name = 'level set variable',
   ...                        mesh = mesh,
   ...                        value = (1,1,1,1,1,-1,-1,-1,-1,-1))
   >>> var.calcDistanceFunction()
   
The problem can then be solved by executing the `solve()` method of the equation.

   >>> if __name__ == '__main__':
   ...     import fipy.viewers
   ...     viewer = fipy.viewers.make(vars = var,
   ...                                limits = {'datamin': -5., 'datamax': 5.})
   ...     viewer.plot()

The result can be tested with the following commands.

   >>> import Numeric
   >>> Numeric.allclose(var, dx * (-Numeric.arange(nx) + (nx - 1) / 2. ))
   1

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    raw_input("finished")
