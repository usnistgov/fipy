#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 3/8/05 {5:05:42 PM} 
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

To run this example from the base FiPy directory, type::
    
    $ examples/diffusion/steadyState/mesh1D/input.py
    
at the command line.  A display of the result should appear and the word 
`finished` in the terminal.

This example takes the user through assembling a simple
problem with FiPy.  It describes a steady 1D diffusion problem with
fixed value boundary conditions such that,

.. raw:: latex

   $$ \nabla \cdot (D \nabla \phi) = 0 $$

with initial conditions

.. raw:: latex

   $\phi = 0$ at $t = 0$,

boundary conditions

.. raw:: latex

   $$ \phi =
   \begin{cases}
       0& \text{at $x = 0$,} \\
       1& \text{at $x = 1$,}
   \end{cases} $$

and parameter value

.. raw:: latex

   $D = 1$.

The first step is to create a mesh with 50 elements. The `Grid1D`
object represents a linear structured grid. The parameter `dx`
refers to the grid spacing (set to unity here).

    >>> nx = 50
    >>> dx = 1.
    >>> from fipy.meshes.numMesh.grid1D import Grid1D
    >>> mesh = Grid1D(nx = nx, dx = dx)

The solution of all equations in FiPy requires a variable. These variables store
values on various parts of the mesh. In this case we need a
`CellVariable` object as the solution is sought on the cell
centers. The boundary conditions are given by `valueLeft = 0` and
`valueRight = 1`. The initial value for the variable is set to `value = valueLeft`.

    >>> valueLeft = 0
    >>> valueRight = 1
    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(name = "solution variable", mesh = mesh, value = valueLeft)

Boundary conditions are given to the equation via a `Tuple`
(list). Boundary conditions are formed with a value and a set of faces
over which they apply. For example here the exterior faces on the left
of the domain are extracted by `mesh.getFacesLeft()`. These faces and
a value (`valueLeft`) are passed to a `FixedValue` boundary
condition. Note that the `FixedFlux(someFaces, 0.)` is the default
boundary condition if no boundary conditions are specified for
exterior faces.

    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> boundaryConditions = (FixedValue(mesh.getFacesRight(),valueRight),
    ...                       FixedValue(mesh.getFacesLeft(),valueLeft))


The steady-state diffusion equation

.. raw:: latex

   $$ \nabla \cdot (D \nabla \phi) = 0 $$
   is represented in \FiPy{} by an `ImplicitDiffusionTerm` object.

    >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    >>> ImplicitDiffusionTerm().solve(var = var, boundaryConditions = boundaryConditions)
    
To test the solution, the analytical result is required. The `x`
coordinates from the mesh are gathered and the length of the domain
`Lx` is calculated.  An array, `analyticalArray`, is calculated to
compare with the numerical result,

    >>> x = mesh.getCellCenters()[:,0]
    >>> Lx = nx * dx
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx

Finally the analytical and numerical results are compared with a
tolerance of `1e-10`.

    >>> var.allclose(analyticalArray, rtol = 1e-10, atol = 1e-10)
    1

The function 'fipy.viewers.make()' returns a suitable viewer depending
on available viewers and the dimension of the mesh.

    >>> if __name__ == '__main__':
    ...     import fipy.viewers
    ...     viewer = fipy.viewers.make(vars = var, limits = {'datamin': 0., 'datamax': 1.})
    ...     viewer.plot()
    ...     raw_input("press key to continue")

..

------

If this example had been written primarily as a script, instead of as
documentation, we would delete every line that does not begin with
either "``>>>``" or "``...``", and then delete those prefixes from the
remaining lines, leaving::
    
    nx = 50
    dx = 1.
    from fipy.meshes.grid2D import Grid2D
    mesh = Grid2D(nx = nx, dx = dx)
    
    valueLeft = 0
    valueRight = 1
    from fipy.variables.cellVariable import CellVariable
    var = CellVariable(name = "solution variable", mesh = mesh, value = valueLeft)
    
    from fipy.boundaryConditions.fixedValue import FixedValue
    boundaryConditions = (FixedValue(mesh.getFacesRight(),valueRight),
			  FixedValue(mesh.getFacesLeft(),valueLeft))

    from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    ImplicitDiffusionTerm().solve(var = var, boundaryConditions = boundaryConditions)
    
    x = mesh.getCellCenters()[:,0]
    Lx = nx * dx
    analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    
    var.allclose(analyticalArray, rtol = 1e-10, atol = 1e-10)
    
    if __name__ == '__main__':
        import fipy.viewers
        viewer = fipy.viewers.make(vars = var, limits = {'datamin': 0., 'datamax': 1.})
        viewer.plot()
	
Your own scripts will tend to look like this, although you can always write
them as doctest scripts if you choose.  You can obtain a plain script
like this from one of the examples by typing::
    
    $ python setup.py copy_script --From examples/.../input.py --To myInput.py

at the command line.

Most of the FiPy examples will be a
mixture of plain scripts and doctest documentation/tests.  
"""

__docformat__ = 'restructuredtext'

def script():
    """
    Return the documentation for this module as a script that can be
    invoked to initialize other scripts.
    """
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.getScript(__name__)

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus.getScript())
    raw_input("finished")
