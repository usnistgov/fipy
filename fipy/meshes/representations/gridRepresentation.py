#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "gridRepresentation.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
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
 # ###################################################################
 ##

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
        >>> fp.dump.write(v, filename='dump.gz')
        >>> v0 = fp.dump.read(filename='dump.gz')
        >>> print (v == v0.mesh.x).all()
        True
        
        >>> m = fp.PeriodicGrid1D(nx=10)
        >>> v = fp.CellVariable(mesh=m, value=m.x)
        >>> fp.dump.write(v, filename='dump.gz')
        >>> v0 = fp.dump.read(filename='dump.gz')
        >>> print (v == v0.mesh.x).all()
        True
        
        >>> m = fp.Tri2D(nx=10, ny=10)
        >>> v = fp.CellVariable(mesh=m, value=m.x)
        >>> fp.dump.write(v, filename='dump.gz')
        >>> v0 = fp.dump.read(filename='dump.gz')
        >>> print (v == v0.mesh.x).all()
        True
        
        >>> m = fp.SkewedGrid2D(nx=10, ny=10)
        >>> v = fp.CellVariable(mesh=m, value=m.x)
        >>> fp.dump.write(v, filename='dump.gz')
        >>> v0 = fp.dump.read(filename='dump.gz')
        >>> print (v == v0.mesh.x).all()
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
