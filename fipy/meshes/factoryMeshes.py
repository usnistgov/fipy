#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Alexander Mont <alexander.mont@nist.gov>
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import parallelComm
from fipy.tools import numerix

__all__ = ["Grid3D", "Grid2D", "Grid1D", "CylindricalGrid2D", "CylindricalGrid1D"]

def _dnl(dx, nx, Lx):
    """
    Initialize arguments for grid classes based on an over determined
    set of initial arguments. The order of precedence is `nx` then
    `Lx` then `dx`. i.e. If `Lx` is specfied the length of the domain
    is always `Lx` regardless of `dx`.

    :Parameters:

      - `dx`: grid spacing
      - `nx`: number of cells
      - `Lx`: the domain length 

    >>> print _dnl(None, None, None)
    (None, 1)
    >>> print _dnl(None, 1.1, 2.)
    (2.0, 1)
    >>> print _dnl(3., None, 7.5)
    (3.75, 2)
    >>> print "(%1.1f, %i)" % _dnl(2.2, 4, None)
    (2.2, 4)
    >>> print _dnl(1., 6, 15.)
    (2.5, 6)
    """
    if Lx is None:
        if nx is None:
            nx = 1
    else:
        if nx is None:
            nx = int(Lx / dx) or 1
        dx = Lx / int(nx)

    return dx, int(nx)

def Grid3D(dx=1., dy=1., dz=1.,
           nx=None, ny=None, nz=None,
           Lx=None, Ly=None, Lz=None,
           overlap=2, communicator=parallelComm):
    
    r""" Factory function to select between UniformGrid3D and
    NonUniformGrid3D.  If `Lx` is specified the length of the domain
    is always `Lx` regardless of `dx`.

    :Parameters:

      - `dx`: grid spacing in the horizontal direction
      - `dy`: grid spacing in the vertical direction
      - `dz`: grid spacing in the z-direction
      - `nx`: number of cells in the horizontal direction
      - `ny`: number of cells in the vertical direction
      - `nz`: number of cells in the z-direction
      - `Lx`: the domain length in the horizontal direction
      - `Ly`: the domain length in the vertical direction
      - `Lz`: the domain length in the z-direction
      - `overlap`: the number of overlapping cells for parallel
        simulations. Generally 2 is adequate. Higher order equations or
        discretizations require more.
      - `communicator`: either `fipy.tools.parallelComm` or
        `fipy.tools.serialComm`. Select `fipy.tools.serialComm` to create a
        serial mesh when running in parallel. Mostly used for test
        purposes.
    
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
    r""" Factory function to select between UniformGrid2D and
    NonUniformGrid2D.  If `Lx` is specified the length of the domain
    is always `Lx` regardless of `dx`.

    :Parameters:

        - `dx`: grid spacing in the horizontal direction
        - `dy`: grid spacing in the vertical direction
        - `nx`: number of cells in the horizontal direction
        - `ny`: number of cells in the vertical direction
        - `Lx`: the domain length in the horizontal direction
        - `Ly`: the domain length in the vertical direction
        - `overlap`: the number of overlapping cells for parallel
          simulations. Generally 2 is adequate. Higher order equations or
          discretizations require more.
        - `communicator`: either `fipy.tools.parallelComm` or
          `fipy.tools.serialComm`. Select `fipy.tools.serialComm` to create a
          serial mesh when running in parallel. Mostly used for test
          purposes.
    
    >>> print Grid2D(Lx=3., nx=2).dx
    1.5

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
    r""" Factory function to select between UniformGrid1D and
    NonUniformGrid1D.  If `Lx` is specified the length of the domain
    is always `Lx` regardless of `dx`.

    :Parameters:

      - `dx`: grid spacing in the horizonal direction
      - `nx`: number of cells in the horizonal direction
      - `Lx`: the domain length in the horizonal direction
      - `overlap`: the number of overlapping cells for parallel
        simulations. Generally 2 is adequate. Higher order equations or
        discretizations require more.
      - `communicator`: either `fipy.tools.parallelComm` or
        `fipy.tools.serialComm`. Select `fipy.tools.serialComm` to create a
        serial mesh when running in parallel. Mostly used for test
        purposes.
    
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
                      origin=((0,),(0,)),
                      overlap=2,
                      communicator=parallelComm):

    r""" Factory function to select between CylindricalUniformGrid2D and
    CylindricalNonUniformGrid2D. If `Lx` is specified the length of
    the domain is always `Lx` regardless of `dx`.

    :Parameters:

      - `dr` or `dx`: grid spacing in the radial direction
      - `dz` or `dy`: grid spacing in the vertical direction
      - `nr` or `nx`: number of cells in the radial direction
      - `nz` or `ny`: number of cells in the vertical direction
      - `Lr` or `Lx`: the domain length in the radial direction
      - `Lz` or `Ly`: the domain length in the vertical direction
      - `origin` : position of the mesh's origin in the form ((x,),(y,))
      - `overlap`: the number of overlapping cells for parallel
        simulations. Generally 2 is adequate. Higher order equations or
        discretizations require more.
      - `communicator`: either `fipy.tools.parallelComm` or
        `fipy.tools.serialComm`. Select `fipy.tools.serialComm` to create a
        serial mesh when running in parallel. Mostly used for test
        purposes.
    
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

    r""" Factory function to select between CylindricalUniformGrid1D and
    CylindricalNonUniformGrid1D. If `Lx` is specified the length of
    the domain is always `Lx` regardless of `dx`.

    :Parameters:

      - `dr` or `dx`: grid spacing in the radial direction
      - `nr` or `nx`: number of cells in the radial direction
      - `Lr` or `Lx`: the domain length in the radial direction
      - `origin` : position of the mesh's origin in the form (x,)
      - `overlap`: the number of overlapping cells for parallel
        simulations. Generally 2 is adequate. Higher order equations or
        discretizations require more.
      - `communicator`: either `fipy.tools.parallelComm` or
        `fipy.tools.serialComm`. Select `fipy.tools.serialComm` to create a
        serial mesh when running in parallel. Mostly used for test
        purposes.
    
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

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
