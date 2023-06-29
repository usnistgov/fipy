r"""

To run this example from the base FiPy directory type::

    $ python examples/diffusion/steadyState/mesh1D/tri2Dinput.py

at the command line. A contour plot should appear and the word `finished`
in the terminal.

This example is similar to the example found in
:mod:`examples.diffusion.mesh1D`, however, the `mesh` is a
:class:`fipy.meshes.tri2D.Tri2D` object rather than a
:func:`~fipy.meshes.factoryMeshes.Grid1D` object.

Here, one time step is executed to implicitly find the steady state
solution.

    >>> DiffusionTerm().solve(var)

To test the solution, the analytical result is required. The `x`
coordinates from the mesh are gathered and the length of the domain,
`Lx`, is calculated.  An array, `analyticalArray`, is calculated to
compare with the numerical result,

    >>> x = mesh.cellCenters[0]
    >>> Lx = nx * dx
    >>> analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx

Finally the analytical and numerical results are compared with a
tolerance of `1e-10`.

    >>> print(var.allclose(analyticalArray))
    1

"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input
from fipy import Tri2D, CellVariable, DiffusionTerm, Viewer

nx = 50
dx = 1.

mesh = Tri2D(dx = dx, nx = nx)

valueLeft = 0.
valueRight = 1.
var = CellVariable(name = "solution-variable", mesh = mesh, value = valueLeft)

var.constrain(valueLeft, mesh.facesLeft)
var.constrain(valueRight, mesh.facesRight)

if __name__ == '__main__':
    DiffusionTerm().solve(var)
    viewer = Viewer(vars=var)
    viewer.plot()
    x = mesh.cellCenters[0]
    Lx = nx * dx
    analyticalArray = valueLeft + (valueRight - valueLeft) * x / Lx
    print(var.allclose(analyticalArray))
    input("finished")

