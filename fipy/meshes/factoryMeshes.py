from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import parallelComm
from fipy.tools import numerix

__all__ = ["Grid3D", "Grid2D", "Grid1D", "CylindricalGrid2D", "CylindricalGrid1D", "SphericalGrid1D"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

def _dnl(dx, nx, Lx):
    """
    Initialize arguments for grid classes based on an over determined set
    of initial arguments.  The order of precedence is `nx` then `Lx` then
    `dx`.  i.e. If `Lx` is specified the length of
    the domain is always `Lx` regardless of `dx`.

    >>> print(_dnl(None, None, None))
    (None, 1)
    >>> print(_dnl(None, 1.1, 2.))
    (2.0, 1)
    >>> print(_dnl(3., None, 7.5))
    (3.75, 2)
    >>> print("(%1.1f, %i)" % _dnl(2.2, 4, None))
    (2.2, 4)
    >>> print(_dnl(1., 6, 15.))
    (2.5, 6)

    Parameters
    ----------
    dx : float
        Grid spacing
    nx : int
        Number of cells
    Lx : float
        The domain length
    """
    if Lx is None:
        if nx is None:
            nx = 1
    else:
        if nx is None:
            nx = Lx // dx or 1
        dx = Lx / int(nx)

    return dx, int(nx)

def Grid3D(dx=1., dy=1., dz=1.,
           nx=None, ny=None, nz=None,
           Lx=None, Ly=None, Lz=None,
           overlap=2, communicator=parallelComm):

    r"""Create a 3D Cartesian mesh

    Factory function to select between
    :class:`~fipy.meshes.uniformGrid3D.UniformGrid3D` and
    :class:`~fipy.meshes.nonUniformGrid3D.NonUniformGrid3D`.  If `L{x,y,z}`
    is specified, the length of the domain is always `L{x,y,z}` regardless
    of `d{x,y,z}`, unless `d{x,y,z}` is a list of spacings, in which case
    `L{x,y,z}` will be the sum of `d{x,y,z}` and `n{x,y,z}` will be the
    count of `d{x,y,z}`.

    Parameters
    ----------
    dx : float
        Grid spacing in the horizontal direction
    dy : float
        Grid spacing in the vertical direction
    dz : float
        Grid spacing in the depth direction
    nx : int
        Number of cells in the horizontal direction
    ny : int
        Number of cells in the vertical direction
    nz : int
        Number of cells in the depth direction
    Lx : float
        Domain length in the horizontal direction
    Ly : float
        Domain length in the vertical direction
    Lz : float
        Domain length in the depth direction
    overlap : int
        Number of overlapping cells for parallel simulations.  Generally 2
        is adequate.  Higher order equations or discretizations require
        more.
    communicator : ~fipy.tools.comms.commWrapper.CommWrapper
        MPI communicator to use.  Select :attr:`~fipy.tools.serialComm` to
        create a serial mesh when running in parallel; mostly used for test
        purposes.  (default: :attr:`~fipy.tools.parallelComm`).
    """

    if numerix.getShape(dx) == () \
      and numerix.getShape(dy) == () \
      and numerix.getShape(dz) == ():

        dx, nx = _dnl(dx, nx, Lx)
        dy, ny = _dnl(dy, ny, Ly)
        dz, nz = _dnl(dz, nz, Lz)
        from fipy.meshes.uniformGrid3D import UniformGrid3D
        return UniformGrid3D(dx = dx, dy = dy, dz = dz,
                             nx = nx or 1, ny = ny or 1, nz = nz or 1,
                             overlap=overlap, communicator=communicator)
    else:
        from fipy.meshes.nonUniformGrid3D import NonUniformGrid3D
        return NonUniformGrid3D(dx = dx, dy = dy, dz = dz, nx = nx, ny = ny, nz = nz,
                                overlap=overlap, communicator=communicator)

def Grid2D(dx=1., dy=1., nx=None, ny=None, Lx=None, Ly=None, overlap=2, communicator=parallelComm):
    r"""Create a 2D Cartesian mesh

    Factory function to select between
    :class:`~fipy.meshes.uniformGrid2D.UniformGrid2D` and
    :class:`~fipy.meshes.nonUniformGrid2D.NonUniformGrid2D`.  If `L{x,y}`
    is specified, the length of the domain is always `L{x,y}` regardless of
    `d{x,y}`, unless `d{x,y}` is a list of spacings, in which case `L{x,y}`
    will be the sum of `d{x,y}` and `n{x,y}` will be the count of `d{x,y}`.

    >>> print(Grid2D(Lx=3., nx=2).dx)
    1.5

    Parameters
    ----------
    dx : float
        Grid spacing in the horizontal direction
    dy : float
        Grid spacing in the vertical direction
    nx : int
        Number of cells in the horizontal direction
    ny : int
        Number of cells in the vertical direction
    Lx : float
        Domain length in the horizontal direction
    Ly : float
        Domain length in the vertical direction
    overlap : int
        Number of overlapping cells for parallel simulations.  Generally 2
        is adequate.  Higher order equations or discretizations require
        more.
    communicator : ~fipy.tools.comms.commWrapper.CommWrapper
        MPI communicator to use.  Select :attr:`~fipy.tools.serialComm` to
        create a serial mesh when running in parallel; mostly used for test
        purposes.  (default: :attr:`~fipy.tools.parallelComm`).
    """

    if numerix.getShape(dx) == () and numerix.getShape(dy) == ():

        dx, nx = _dnl(dx, nx, Lx)
        dy, ny = _dnl(dy, ny, Ly)

        from fipy.meshes.uniformGrid2D import UniformGrid2D
        return UniformGrid2D(dx=dx, dy=dy,
                             nx=nx, ny=ny,
                             overlap=overlap,
                             communicator=communicator)
    else:
        from fipy.meshes.nonUniformGrid2D import NonUniformGrid2D
        return NonUniformGrid2D(dx=dx, dy=dy, nx=nx, ny=ny, overlap=overlap, communicator=communicator)

def Grid1D(dx=1., nx=None, Lx=None, overlap=2, communicator=parallelComm):
    r"""Create a 1D Cartesian mesh

    Factory function to select between
    :class:`~fipy.meshes.uniformGrid1D.UniformGrid1D` and
    :class:`~fipy.meshes.nonUniformGrid1D.NonUniformGrid1D`.  If `Lx` is
    specified the length of the domain is always `Lx` regardless of `dx`,
    unless `dx` is a list of spacings, in which case `Lx` will be the sum
    of `dx` and `nx` will be the count of `dx`.

    Parameters
    ----------
    dx : float
        Grid spacing in the horizontal direction
    nx : int
        Number of cells in the horizontal direction
    Lx : float
        Domain length in the horizontal direction
    overlap : int
        Number of overlapping cells for parallel simulations.  Generally 2
        is adequate.  Higher order equations or discretizations require
        more.
    communicator : ~fipy.tools.comms.commWrapper.CommWrapper
        MPI communicator to use.  Select :attr:`~fipy.tools.serialComm` to
        create a serial mesh when running in parallel; mostly used for test
        purposes.  (default: :attr:`~fipy.tools.parallelComm`).
    """

    if numerix.getShape(dx) == ():
        dx, nx = _dnl(dx, nx, Lx)
        from fipy.meshes.uniformGrid1D import UniformGrid1D
        return UniformGrid1D(dx=dx, nx=nx, overlap=overlap, communicator=communicator)
    else:
        from fipy.meshes.nonUniformGrid1D import NonUniformGrid1D
        return NonUniformGrid1D(dx=dx, nx=nx, overlap=overlap, communicator=communicator)

def CylindricalGrid2D(dr=None, dz=None,
                      nr=None, nz=None,
                      Lr=None, Lz=None,
                      dx=1., dy=1.,
                      nx=None, ny=None,
                      Lx=None, Ly=None,
                      origin=((0,), (0,)),
                      overlap=2,
                      communicator=parallelComm):

    r"""Create a 2D cylindrical mesh

    Factory function to select between
    :class:`~fipy.meshes.cylindricalUniformGrid2D.CylindricalUniformGrid2D`
    and
    :class:`~fipy.meshes.cylindricalNonUniformGrid2D.CylindricalNonUniformGrid2D`.
    If `Lr` is specified the length of the domain is always `Lr` regardless
    of `dr`, unless `dr` is a list of spacings, in which case `Lr` will be
    the sum of `dr`.

    Parameters
    ----------
    dr : float
        Grid spacing in the radial direction. Alternative: `dx`.
    dz : float
        grid spacing in the vertical direction. Alternative: `dy`.
    nr : int
        Number of cells in the radial direction. Alternative: `nx`.
    nz : int
        Number of cells in the vertical direction. Alternative: `ny`.
    Lr : float
        Domain length in the radial direction. Alternative: `Lx`.
    Lz : float
        Domain length in the vertical direction. Alternative: `Ly`.
    overlap : int
        the number of overlapping cells for parallel
        simulations. Generally 2 is adequate. Higher order equations or
        discretizations require more.
    communicator : ~fipy.tools.comms.commWrapper.CommWrapper
        MPI communicator to use.  Select :attr:`~fipy.tools.serialComm` to
        create a serial mesh when running in parallel; mostly used for test
        purposes.  (default: :attr:`~fipy.tools.parallelComm`).
    """
    if dr is not None:
        dx = dr

    if dz is not None:
        dy = dz


    nx = nr or nx
    ny = nz or ny

    Lx = Lr or Lx
    Ly = Lz or Ly

    if numerix.getShape(dx) == () and numerix.getShape(dy) == ():

        dx, nx = _dnl(dx, nx, Lx)
        dy, ny = _dnl(dy, ny, Ly)
        from fipy.meshes.cylindricalUniformGrid2D import CylindricalUniformGrid2D
        return CylindricalUniformGrid2D(dx=dx, dy=dy, nx=nx or 1, ny=ny or 1, origin=origin, overlap=overlap, communicator=communicator)
    else:
        from fipy.meshes.cylindricalNonUniformGrid2D import CylindricalNonUniformGrid2D
        return CylindricalNonUniformGrid2D(dx=dx, dy=dy, nx=nx, ny=ny, origin=origin, overlap=overlap, communicator=communicator)

def CylindricalGrid1D(dr=None, nr=None, Lr=None,
                      dx=1., nx=None, Lx=None,
                      origin=(0,), overlap=2, communicator=parallelComm):

    r"""Create a 2D cylindrical mesh

    Factory function to select between
    :class:`~fipy.meshes.cylindricalUniformGrid1D.CylindricalUniformGrid1D`
    and
    :class:`~fipy.meshes.cylindricalNonUniformGrid1D.CylindricalNonUniformGrid1D`.
    If `Lr` is specified the length of the domain is always `Lr` regardless
    of `dr`, unless `dr` is a list of spacings, in which case `Lr` will be
    the sum of `dr`.

    Parameters
    ----------
    dr : float
        Grid spacing in the radial direction. Alternative: `dx`.
    nr : int
        Number of cells in the radial direction. Alternative: `nx`.
    Lr : float
        Domain length in the radial direction. Alternative: `Lx`.
    overlap : int
        the number of overlapping cells for parallel
        simulations. Generally 2 is adequate. Higher order equations or
        discretizations require more.
    communicator : ~fipy.tools.comms.commWrapper.CommWrapper
        MPI communicator to use.  Select :attr:`~fipy.tools.serialComm` to
        create a serial mesh when running in parallel; mostly used for test
        purposes.  (default: :attr:`~fipy.tools.parallelComm`).
    """

    if dr is not None:
        dx = dr

    nx = nr or nx
    Lx = Lr or Lx

    if numerix.getShape(dx) == ():
        dx, nx = _dnl(dx, nx, Lx)
        from fipy.meshes.cylindricalUniformGrid1D import CylindricalUniformGrid1D
        return CylindricalUniformGrid1D(dx=dx, nx=nx or 1, origin=origin, overlap=overlap, communicator=parallelComm)
    else:
        from fipy.meshes.cylindricalNonUniformGrid1D import CylindricalNonUniformGrid1D
        return CylindricalNonUniformGrid1D(dx=dx, nx=nx, origin=origin, overlap=overlap, communicator=parallelComm)

def SphericalGrid1D(dr=None, nr=None, Lr=None,
                    dx=1., nx=None, Lx=None,
                    origin=(0,), overlap=2, communicator=parallelComm):

    r"""Create a 1D spherical mesh

    Factory function to select between
    :class:`~fipy.meshes.sphericalUniformGrid1D.SphericalUniformGrid1D` and
    :class:`~fipy.meshes.sphericalNonUniformGrid1D.SphericalNonUniformGrid1D`.
    If `Lr` is specified the length of the domain is always `Lr` regardless
    of `dr`, unless `dr` is a list of spacings, in which case `Lr` will be
    the sum of `dr`.

    Parameters
    ----------
    dr : float
        Grid spacing in the radial direction. Alternative: `dx`.
    nr : int
        Number of cells in the radial direction. Alternative: `nx`.
    Lr : float
        Domain length in the radial direction. Alternative: `Lx`.
    overlap : int
        the number of overlapping cells for parallel
        simulations. Generally 2 is adequate. Higher order equations or
        discretizations require more.
    communicator : ~fipy.tools.comms.commWrapper.CommWrapper
        MPI communicator to use.  Select :attr:`~fipy.tools.serialComm` to
        create a serial mesh when running in parallel; mostly used for test
        purposes.  (default: :attr:`~fipy.tools.parallelComm`).
    """

    if dr is not None:
        dx = dr

    nx = nr or nx
    Lx = Lr or Lx

    if numerix.getShape(dx) == ():
        dx, nx = _dnl(dx, nx, Lx)
        from fipy.meshes.sphericalUniformGrid1D import SphericalUniformGrid1D
        return SphericalUniformGrid1D(dx=dx, nx=nx or 1, origin=origin, overlap=overlap, communicator=parallelComm)
    else:
        from fipy.meshes.sphericalNonUniformGrid1D import SphericalNonUniformGrid1D
        return SphericalNonUniformGrid1D(dx=dx, nx=nx, origin=origin, overlap=overlap, communicator=parallelComm)
    
def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()

