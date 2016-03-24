#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "sphere.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://ctcms.nist.gov
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

"""
An interesting problem is to solve an equation on a 2D geometry that
is embedded in 3D space, such as diffusion on the surface of a sphere
(with nothing either inside or outside the sphere). This example
demonstrates how to create the required mesh.

    >>> from fipy import Gmsh2DIn3DSpace, CellVariable, MayaviClient
    >>> from fipy.tools import numerix

    >>> mesh = Gmsh2DIn3DSpace('''
    ...     radius = 5.0;
    ...     cellSize = 0.3;
    ...
    ...     // create inner 1/8 shell
    ...     Point(1) = {0, 0, 0, cellSize};
    ...     Point(2) = {-radius, 0, 0, cellSize};
    ...     Point(3) = {0, radius, 0, cellSize};
    ...     Point(4) = {0, 0, radius, cellSize};
    ...     Circle(1) = {2, 1, 3};
    ...     Circle(2) = {4, 1, 2};
    ...     Circle(3) = {4, 1, 3};
    ...     Line Loop(1) = {1, -3, 2} ;
    ...     Ruled Surface(1) = {1};
    ...
    ...     // create remaining 7/8 inner shells
    ...     t1[] = Rotate {{0,0,1},{0,0,0},Pi/2} {Duplicata{Surface{1};}};
    ...     t2[] = Rotate {{0,0,1},{0,0,0},Pi} {Duplicata{Surface{1};}};
    ...     t3[] = Rotate {{0,0,1},{0,0,0},Pi*3/2} {Duplicata{Surface{1};}};
    ...     t4[] = Rotate {{0,1,0},{0,0,0},-Pi/2} {Duplicata{Surface{1};}};
    ...     t5[] = Rotate {{0,0,1},{0,0,0},Pi/2} {Duplicata{Surface{t4[0]};}};
    ...     t6[] = Rotate {{0,0,1},{0,0,0},Pi} {Duplicata{Surface{t4[0]};}};
    ...     t7[] = Rotate {{0,0,1},{0,0,0},Pi*3/2} {Duplicata{Surface{t4[0]};}};
    ...
    ...     // create entire inner and outer shell
    ...     Surface Loop(100)={1,t1[0],t2[0],t3[0],t7[0],t4[0],t5[0],t6[0]};
    ... ''').extrude(extrudeFunc=lambda r: 1.1 * r) # doctest: +GMSH

    >>> x, y, z = mesh.cellCenters # doctest: +GMSH

    >>> var = CellVariable(mesh=mesh, value=x * y * z, name="x*y*z") # doctest: +GMSH

    >>> if __name__ == '__main__':
    ...     viewer = MayaviClient(vars=var)
    ...     viewer.plot()

   >>> max(numerix.sqrt(x**2 + y**2 + z**2)) < 5.3 # doctest: +GMSH
   True
   >>> min(numerix.sqrt(x**2 + y**2 + z**2)) > 5.2 # doctest: +GMSH
   True

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')
