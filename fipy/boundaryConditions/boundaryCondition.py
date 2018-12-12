__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.variables.variable import Variable
from fipy.tools.dimensions.physicalField import PhysicalField

__all__ = ["BoundaryCondition"]

class BoundaryCondition(object):
    """
    Generic boundary condition base class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self,faces,value):
        """
        :Parameters:
            - `faces`: A `list` or `tuple` of exterior `Face` objects to which this condition applies.
            - `value`: The value to impose.
        """
        if self.__class__ is BoundaryCondition:
            raise NotImplementedError, "can't instantiate abstract base class"

        self.faces = faces
        if not (isinstance(value, PhysicalField) or isinstance(value, Variable)):
            value = PhysicalField(value)
        self.value = value

        if not (self.faces | self.faces.mesh.exteriorFaces
                == self.faces.mesh.exteriorFaces).value.all():
            raise IndexError, 'Face list has interior faces'

        self.adjacentCellIDs = self.faces.mesh._adjacentCellIDs[0][self.faces.value]
        self.boundaryConditionApplied = False

    def _buildMatrix(self, SparseMatrix, Ncells, MaxFaces, coeff):
        """Return the effect of this boundary condition on the equation
        solution matrices.

        `_buildMatrix()` is called by each `Term` of each `Equation`.

        :Parameters:
          - `SparseMatrix`: Sparse matrix class to use
          - `Ncells`:       Number of cells (to build **L** and **b**)
          - `MaxFaces`:     Maximum number of faces per cell (to build **L**)
          - `coeff`:        Contribution due to this face

        A `tuple` of (`LL`, `bb`) is calculated, to be added to the Term's
        (**L**, **b**) matrices.
        """
        raise NotImplementedError

    def _getDerivative(self, order):
        """Return a tuple of the boundary conditions to apply
        to the term and to the derivative of the term
        """
        if order == 0:
            return self
        else:
            return None

    def __repr__(self):
        return "%s(faces = %s, value = %s)" % (self.__class__.__name__, `self.faces`, `self.value`)

    def _resetBoundaryConditionApplied(self):
        self.boundaryConditionApplied = False

    def _test(self):
        """
        The `BoundaryCondition` class should raise an error when
        invoked with internal faces. Don't use the `BoundaryCondition`
        class in this manner. This is merely a test.

        >>> from fipy.meshes import Grid1D
        >>> mesh = Grid1D(nx = 2)
        >>> from fipy.tools import parallelComm
        >>> if parallelComm.procID == 0:
        ...     bc = __BoundaryCondition(mesh.interiorFaces, 0)
        ... else:
        ...     raise IndexError("Face list has interior faces")
        Traceback (most recent call last):
            ...
        IndexError: Face list has interior faces
        """
        pass

class __BoundaryCondition(BoundaryCondition):
    """
    Dummy subclass for tests
    """
    pass

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
