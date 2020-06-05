"""
1D Mesh
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import parallelComm

from fipy.meshes.nonUniformGrid1D import NonUniformGrid1D

__all__ = ["SphericalNonUniformGrid1D"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class SphericalNonUniformGrid1D(NonUniformGrid1D):
    """
    Creates a 1D spherical grid mesh.

        >>> mesh = SphericalNonUniformGrid1D(nx = 3)
        >>> print(mesh.cellCenters)
        [[ 0.5  1.5  2.5]]

        >>> mesh = SphericalNonUniformGrid1D(dx = (1, 2, 3))
        >>> print(mesh.cellCenters)
        [[ 0.5  2.   4.5]]

        >>> print(numerix.allclose(mesh.cellVolumes, (0.5, 13., 94.5))) # doctest: +PROCESSOR_0
        True

        >>> mesh = SphericalNonUniformGrid1D(nx = 2, dx = (1, 2, 3))
        Traceback (most recent call last):
        ...
        IndexError: nx != len(dx)

        >>> mesh = SphericalNonUniformGrid1D(nx=2, dx=(1., 2.)) + ((1.,),)
        >>> print(mesh.cellCenters)
        [[ 1.5  3. ]]
        >>> print(numerix.allclose(mesh.cellVolumes, (3.5, 28))) # doctest: +PROCESSOR_0
        True

    """
    def __init__(self, dx=1., nx=None, origin=(0,), overlap=2, communicator=parallelComm, *args, **kwargs):
        scale = PhysicalField(value=1, unit=PhysicalField(value=dx).unit)
        self.origin = PhysicalField(value=origin)
        self.origin /= scale

        super(SphericalNonUniformGrid1D, self).__init__(dx=dx,
                                                        nx=nx,
                                                        overlap=overlap,
                                                        communicator=communicator,
                                                        *args,
                                                        **kwargs)

        self.vertexCoords += origin
        self.args['origin'] = origin

    def _calcFaceCenters(self):
        faceCenters = super(SphericalNonUniformGrid1D, self)._calcFaceCenters()
        return faceCenters + self.origin

    def _calcFaceAreas(self):
        return self._calcFaceCenters()[0] * self._calcFaceCenters()[0]

    def _calcCellVolumes(self):
        return super(SphericalNonUniformGrid1D, self)._calcCellVolumes() / 2.

    def _translate(self, vector):
        return SphericalNonUniformGrid1D(dx=self.args['dx'], nx=self.args['nx'],
                                         origin=numerix.array(self.args['origin']) + vector,
                                         overlap=self.args['overlap'])

    def __mul__(self, factor):
        return SphericalNonUniformGrid1D(dx=self.args['dx'] * factor, nx=self.args['nx'],
                                         origin=numerix.array(self.args['origin']) * factor,
                                         overlap=self.args['overlap'])

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected. Fixed a bug where the following throws
        an error on solve() when `nx` is a float.

            >>> from fipy import CellVariable, DiffusionTerm
            >>> mesh = SphericalNonUniformGrid1D(nx=3., dx=(1., 2., 3.))
            >>> var = CellVariable(mesh=mesh)
            >>> var.constrain(0., where=mesh.facesRight)
            >>> DiffusionTerm().solve(var)

        This test is for https://github.com/usnistgov/fipy/issues/372. Cell
        volumes were being returned as `binOps` rather than arrays.

            >>> m = SphericalNonUniformGrid1D(dx=(1., 2., 3., 4.), nx=4)
            >>> print(isinstance(m.cellVolumes, numerix.ndarray))
            True
            >>> print(isinstance(m._faceAreas, numerix.ndarray))
            True

        If the above types aren't correct, the divergence operator's value can be a `binOp`

            >>> print(isinstance(CellVariable(mesh=m).arithmeticFaceValue.divergence.value, numerix.ndarray))
            True

        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

