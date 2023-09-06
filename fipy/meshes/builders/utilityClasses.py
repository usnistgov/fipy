from __future__ import unicode_literals
from builtins import object
from builtins import range
from builtins import zip
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import numerix

class _UniformOrigin(object):
    """
    Used to calculate the origin for uniform grids in a
    dimensionally-independent way.
    """

    @staticmethod
    def calcOrigin(origin, offset, ds, scale):
        newOrigin  = PhysicalField(value=origin)
        newOrigin /= scale

        if type(offset) in [int, float]:
            newOrigin += offset * ds[0]
        else:
            newOrigin += [[o*float(d)] for o, d in zip(offset, ds)]

        return newOrigin

class _DOffsets(object):
    """
    For use by non-uniform grid builders.
    """

    @staticmethod
    def calcDOffsets(ds, ns, offset):
        """
        Parameters
        ----------
        ds : list
            Spacing in each grid direction, e.g. `[dx, dy]`
        ns : list
            Number of grid spacings in each direction, e.g. `[nx, ny]`
        offset : list
            Displacement of grid spacings, e.g., `[Ox, Oy]`

        Returns
        -------
        offsetList : list
            Contains the analogous scalars to `XOffset`, `YOffset`, and
            `ZOffset`, whichever are applicable for the dimensionality.
        newDs : list
            proper `[dx, [dy, ...]]` values
        """
        offsetList = []
        newDs = []

        if type(offset) in [int, float]:
            offset = [offset]

        for d, n, i in zip(ds, ns, list(range(len(ds)))):
            if len(numerix.getShape(d)) > 0:
                offsetList.append(numerix.sum(d[0:offset[i]]))
                newDs.append(d[offset[i]:offset[i] + n])
            else:
                if len(offset) == 1:
                    offsetList.append(d * offset[0])
                else:
                    offsetList.append(d * offset[i])

                newDs.append(d)

        return offsetList, newDs

class _AbstractNumPts(object):
    """
    Interface definition for `NumPts` calculators.
    """

    @staticmethod
    def calcNs(ns, ds):
        raise NotImplementedError

class _NonuniformNumPts(_AbstractNumPts):
    """
    For use by non-uniform grid builders.
    """

    @staticmethod
    def calcNs(ns, ds):
        axis = ["x", "y", "z"][:len(ns)]
        newNs = []

        for a, d, n in zip(axis, ds, ns):
            newNs.append(_NonuniformNumPts._calcNumPts(d=d, n=n, axis=a))

        return newNs

    @staticmethod
    def _calcNumPts(d, n = None, axis = "x"):
        """
        Calculate the number of cells along the specified axis, based
        on either the specified number or on the number elements in the
        cell  `d` spacings.

        Used by the `Grid` meshes.

        This tests a bug that was occurring with `PeriodicGrid1D` when
        using a numpy float as the argument for the grid spacing.

           >>> from fipy.meshes.periodicGrid1D import PeriodicGrid1D
           >>> PeriodicGrid1D(nx=2, dx=numerix.float32(1.))
           PeriodicGrid1D(dx=1.0, nx=2)
        """

        if type(d) in [int, float] or numerix.shape(d) == ():
            n = int(n or 1)
        else:
            n = int(n or len(d))
            if n != len(d) and len(d) != 1:
                raise IndexError("n%s != len(d%s)" % (axis, axis))

        return n

class _UniformNumPts(_AbstractNumPts):
    """
    For use by uniform grid builders.
    """

    @staticmethod
    def calcNs(ns, ds):
        return [int(x) for x in ns]

if __name__ == '__main__':
    import doctest
    doctest.testmod()
