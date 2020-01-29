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
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')
