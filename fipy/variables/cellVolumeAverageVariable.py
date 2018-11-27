


__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix

from fipy.variables.variable import Variable
from fipy.variables.cellVariable import CellVariable

class _CellVolumeAverageVariable(Variable):
    """
    Takes a `CellVariable` and evaluates its volume average over all the
    cells.

        >>> from fipy.meshes import Grid2D
        >>> mesh = Grid2D(nx = 2, ny = 2, dx = 2., dy = 5.)
        >>> from fipy.variables.cellVariable import CellVariable
        >>> var = CellVariable(value = (1, 2, 3 ,4), mesh = mesh)
        >>> print _CellVolumeAverageVariable(var)
        2.5

    """
    def __init__(self, var):
        Variable.__init__(self, unit = var.unit)
        self.var = self._requires(var)

    def _calcValue(self):
        mesh = self.var.mesh
        volumes = CellVariable(mesh=mesh, value=mesh.cellVolumes)
        return (self.var * volumes).sum() / volumes.sum()

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
