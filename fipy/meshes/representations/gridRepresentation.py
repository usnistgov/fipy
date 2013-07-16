from __future__ import unicode_literals
from builtins import str
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.meshes.representations.abstractRepresentation import _AbstractRepresentation

class _GridRepresentation(_AbstractRepresentation):

    def getstate(self):
        """Collect the necessary information to ``pickle`` the `Grid` to persistent storage.
        """
        args = self.mesh.args.copy()
        args["_RepresentationClass"] = self.__class__
        return args

    @staticmethod
    def setstate(mesh, state):
        """Create a new `Grid` from ``pickled`` persistent storage.
        """
        mesh.__init__(**state)

    def _repr(self, dns):
        dnstr = []
        for d, n in dns:
            dnstr.append(d + "=" + str(self.mesh.args[d]))
            if self.mesh.args[n] is not None:
                dnstr.append(n + "=" + str(self.mesh.args[n]))

        return "%s(%s)" % (self.mesh.__class__.__name__, ", ".join(dnstr))

    def _test(self):
        """

        Check that the following grid classes can be pickled and unpickled.

        >>> import fipy as fp

        >>> m = fp.PeriodicGrid2DLeftRight(nx=10, ny=10)
        >>> v = fp.CellVariable(mesh=m, value=m.x)
        >>> (f, filename) = fp.dump.write(v, extension='.gz')
        >>> v0 = fp.dump.read(filename, f)
        >>> print((v == v0.mesh.x).all())
        True

        >>> m = fp.PeriodicGrid1D(nx=10)
        >>> v = fp.CellVariable(mesh=m, value=m.x)
        >>> (f, filename) = fp.dump.write(v, extension='.gz')
        >>> v0 = fp.dump.read(filename, f)
        >>> print((v == v0.mesh.x).all())
        True

        >>> m = fp.Tri2D(nx=10, ny=10)
        >>> v = fp.CellVariable(mesh=m, value=m.x)
        >>> (f, filename) = fp.dump.write(v, extension='.gz')
        >>> v0 = fp.dump.read(filename, f)
        >>> print((v == v0.mesh.x).all())
        True

        >>> m = fp.SkewedGrid2D(nx=10, ny=10)
        >>> v = fp.CellVariable(mesh=m, value=m.x)
        >>> (f, filename) = fp.dump.write(v, extension='.gz')
        >>> v0 = fp.dump.read(filename, f)
        >>> print((v == v0.mesh.x).all())
        True

        """

class _Grid1DRepresentation(_GridRepresentation):

    def repr(self):
        return self._repr(dns=[("dx", "nx")])

class _Grid2DRepresentation(_GridRepresentation):

    def repr(self):
        return self._repr(dns=[("dx", "nx"), ("dy", "ny")])

class _Grid3DRepresentation(_GridRepresentation):

    def repr(self):
        return self._repr(dns=[("dx", "nx"), ("dy", "ny"), ("dz", "nz")])

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

