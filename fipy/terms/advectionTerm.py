from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = ["AdvectionTerm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

from fipy.tools.numerix import MA
from fipy.tools import numerix

from fipy.terms.firstOrderAdvectionTerm import FirstOrderAdvectionTerm

class AdvectionTerm(FirstOrderAdvectionTerm):
    r"""

    The `AdvectionTerm` object constructs the `b` vector contribution for
    the advection term given by

    .. math::

       u \abs{\nabla \phi}

    from the advection equation given by:

    .. math::

       \frac{\partial \phi}{\partial t} + u \abs{\nabla \phi} = 0

    The construction of the gradient magnitude term requires upwinding as in the standard
    `FirstOrderAdvectionTerm`. The higher order terms are incorporated as follows.
    The formula used here is given by:

    .. math::

       u_P \abs{\nabla \phi}_P = \max \left( u_P , 0 \right) \left[  \sum_A \min \left( D_{AP}, 0 \right)^2 \right]^{1/2} +  \min \left( u_P , 0 \right) \left[  \sum_A \max \left( D_{AP}, 0 \right)^2 \right]^{1/2}

    where,

    .. math::

       D_{AP} = \frac{ \phi_A - \phi_P } { d_{AP}} - \frac{ d_{AP} } {2} m \left(L_A, L_P \right)

    and

    .. math::

       m\left(x, y\right) &= x \qquad \text{if $\abs{x} \le \abs{y} \forall xy \ge 0$} \\
       m\left(x, y\right) &= y \qquad \text{if $\abs{x} > \abs{y} \forall xy \ge 0$} \\
       m\left(x, y\right) &= 0 \qquad \text{if $xy < 0$}

    also,

    .. math::

       L_A &= \frac{\phi_{AA} + \phi_P - 2 \phi_A}{d_{AP}^2} \\
       L_P &= \frac{\phi_{A} + \phi_{PP} - 2 \phi_P}{d_{AP}^2}

    Here are some simple test cases for this problem:

    >>> from fipy.meshes import Grid1D
    >>> from fipy.solvers import *
    >>> SparseMatrix = LinearPCGSolver()._matrixClass
    >>> mesh = Grid1D(dx = 1., nx = 3)

    Trivial test:

    >>> from fipy.variables.cellVariable import CellVariable
    >>> coeff = CellVariable(mesh = mesh, value = numerix.zeros(3, 'd'))
    >>> v, L, b = AdvectionTerm(0.)._buildMatrix(coeff, SparseMatrix)
    >>> print(numerix.allclose(b, numerix.zeros(3, 'd'), atol = 1e-10)) # doctest: +PROCESSOR_0
    True

    Less trivial test:

    >>> coeff = CellVariable(mesh = mesh, value = numerix.arange(3))
    >>> v, L, b = AdvectionTerm(1.)._buildMatrix(coeff, SparseMatrix)
    >>> print(numerix.allclose(b, numerix.array((0., -1., -1.)), atol = 1e-10)) # doctest: +PROCESSOR_0
    True

    Even less trivial

    >>> coeff = CellVariable(mesh = mesh, value = numerix.arange(3))
    >>> v, L, b = AdvectionTerm(-1.)._buildMatrix(coeff, SparseMatrix)
    >>> print(numerix.allclose(b, numerix.array((1., 1., 0.)), atol = 1e-10)) # doctest: +PROCESSOR_0
    True

    Another trivial test case (more trivial than a trivial test case
    standing on a harpsichord singing "trivial test cases are here again")

    >>> vel = numerix.array((-1, 2, -3))
    >>> coeff = CellVariable(mesh = mesh, value = numerix.array((4, 6, 1)))
    >>> v, L, b = AdvectionTerm(vel)._buildMatrix(coeff, SparseMatrix)
    >>> print(numerix.allclose(b, -vel * numerix.array((2, numerix.sqrt(5**2 + 2**2), 5)), atol = 1e-10)) # doctest: +PROCESSOR_0
    True

    Somewhat less trivial test case:

    >>> from fipy.meshes import Grid2D
    >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 2)
    >>> vel = numerix.array((3, -5, -6, -3))
    >>> coeff = CellVariable(mesh = mesh, value = numerix.array((3, 1, 6, 7)))
    >>> v, L, b = AdvectionTerm(vel)._buildMatrix(coeff, SparseMatrix)
    >>> answer = -vel * numerix.array((2, numerix.sqrt(2**2 + 6**2), 1, 0))
    >>> print(numerix.allclose(b, answer, atol = 1e-10)) # doctest: +PROCESSOR_0
    True

    For the above test cases the `AdvectionTerm` gives the
    same result as the `AdvectionTerm`. The following test imposes a quadratic
    field. The higher order term can resolve this field correctly.

    .. math::

       \phi = x^2

    The returned vector ``b`` should have the value:

    .. math::

       -\abs{\nabla \phi} = -\left|\frac{\partial \phi}{\partial x}\right| = - 2 \abs{x}

    Build the test case in the following way,

    >>> mesh = Grid1D(dx = 1., nx = 5)
    >>> vel = 1.
    >>> coeff = CellVariable(mesh = mesh, value = mesh.cellCenters[0]**2)
    >>> v, L, b = __AdvectionTerm(vel)._buildMatrix(coeff, SparseMatrix)

    The first order term is not accurate. The first and last element are ignored because they
    don't have any neighbors for higher order evaluation

    >>> print(numerix.allclose(CellVariable(mesh=mesh,
    ... value=b).globalValue[1:-1], -2 * mesh.cellCenters.globalValue[0][1:-1]))
    False

    The higher order term is spot on.

    >>> v, L, b = AdvectionTerm(vel)._buildMatrix(coeff, SparseMatrix)
    >>> print(numerix.allclose(CellVariable(mesh=mesh,
    ... value=b).globalValue[1:-1], -2 * mesh.cellCenters.globalValue[0][1:-1]))
    True

    The `AdvectionTerm` will also resolve a circular field with
    more accuracy,

    .. math::

       \phi = \left( x^2 + y^2 \right)^{1/2}

    Build the test case in the following way,

    >>> mesh = Grid2D(dx = 1., dy = 1., nx = 10, ny = 10)
    >>> vel = 1.
    >>> x, y = mesh.cellCenters
    >>> r = numerix.sqrt(x**2 + y**2)
    >>> coeff = CellVariable(mesh = mesh, value = r)
    >>> v, L, b = __AdvectionTerm(1.)._buildMatrix(coeff, SparseMatrix)
    >>> error = CellVariable(mesh=mesh, value=b + 1)
    >>> ans = CellVariable(mesh=mesh, value=b + 1)
    >>> ans[(x > 2) & (x < 8) & (y > 2) & (y < 8)] = 0.123105625618
    >>> print((error <= ans).all())
    True

    The maximum error is large (about 12 %) for the first order advection.

    >>> v, L, b = AdvectionTerm(1.)._buildMatrix(coeff, SparseMatrix)
    >>> error = CellVariable(mesh=mesh, value=b + 1)
    >>> ans = CellVariable(mesh=mesh, value=b + 1)
    >>> ans[(x > 2) & (x < 8) & (y > 2) & (y < 8)] = 0.0201715476598
    >>> print((error <= ans).all())
    True

    The maximum error is 2 % when using a higher order contribution.

    """
    def _getDifferences(self, adjacentValues, cellValues, oldArray, cellToCellIDs, mesh):

        dAP = mesh._cellToCellDistances

##        adjacentGradient = numerix.take(oldArray.grad, cellToCellIDs)
        adjacentGradient = numerix.take(oldArray.grad, mesh._cellToCellIDs, axis=-1)
        adjacentNormalGradient = numerix.dot(adjacentGradient, mesh._cellNormals)
        adjacentUpValues = cellValues + 2 * dAP * adjacentNormalGradient

        cellIDs = numerix.repeat(numerix.arange(mesh.numberOfCells)[numerix.newaxis, ...],
                mesh._maxFacesPerCell, axis=0)
        cellIDs = MA.masked_array(cellIDs, mask = MA.getmask(mesh._cellToCellIDs))
        cellGradient = numerix.take(oldArray.grad, cellIDs, axis=-1)
        cellNormalGradient = numerix.dot(cellGradient, mesh._cellNormals)
        cellUpValues = adjacentValues - 2 * dAP * cellNormalGradient

        cellLaplacian = (cellUpValues + adjacentValues - 2 * cellValues) / dAP**2

        adjacentLaplacian = (adjacentUpValues + cellValues - 2 * adjacentValues) / dAP**2
        adjacentLaplacian = adjacentLaplacian.filled(0)
        cellLaplacian = cellLaplacian.filled(0)

        mm = numerix.where(cellLaplacian * adjacentLaplacian < 0.,
                           0.,
                           numerix.where(abs(cellLaplacian) > abs(adjacentLaplacian),
                                         adjacentLaplacian,
                                         cellLaplacian))

        return FirstOrderAdvectionTerm._getDifferences(self, adjacentValues, cellValues, oldArray, cellToCellIDs, mesh) -  mm * dAP / 2.

class __AdvectionTerm(FirstOrderAdvectionTerm):
    """
    Dummy subclass for tests
    """
    pass

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()


