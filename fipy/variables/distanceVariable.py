#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "distanceVariable.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.variables.cellVariable import CellVariable

from fipy.tests.doctestPlus import register_skipper
import sys
import os

def module_exists(module_name):
    try:
        __import__(module_name)
    except ImportError:
        return False
    else:
        return True

def _checkForLSMLIB():
    return module_exists('pylsmlib')

def _checkForSKFMM():
    return module_exists('skfmm')

def _parseLSMSolver():
    args = [s.lower() for s in sys.argv[1:]]
    # any command-line specified solver takes precedence over environment variables
    if '--lsmlib' in args:
        if _checkForLSMLIB():
            return "lsmlib"
        else:
            return None
    elif '--skfmm' in args:
        if _checkForSKFMM():
            return "skfmm"
        else:
            return None
    elif 'FIPY_LSM' in os.environ:
        return os.environ['FIPY_LSM'].lower()
    elif _checkForLSMLIB():
        return 'lsmlib'
    elif _checkForSKFMM():
        return 'skfmm'
    else:
        return None

LSM_SOLVER = _parseLSMSolver()

register_skipper(flag="LSM",
                 test=lambda : LSM_SOLVER is not None,
                 why="neither `lsmlib` nor `skfmm` can be found on the $PATH")

register_skipper(flag="LSMLIB",
                 test=lambda : LSM_SOLVER == 'lsmlib',
                 why="`lsmlib` must be used to run some tests")

register_skipper(flag="SKFMM",
                 test=lambda : LSM_SOLVER == 'skfmm',
                 why="`skfmm` must be used to run some tests")


__all__ = ["DistanceVariable"]

class DistanceVariable(CellVariable):
    r"""
    A `DistanceVariable` object calculates :math:`\phi` so it satisfies,

    .. math::

       \abs{\nabla \phi} = 1

    using the fast marching method with an initial condition defined
    by the zero level set.  The solution can either be first or second
    order.

    Here we will define a few test cases. Firstly a 1D test case

    >>> from fipy.meshes import Grid1D
    >>> from fipy.tools import serialComm
    >>> mesh = Grid1D(dx = .5, nx = 8, communicator=serialComm)
    >>> from distanceVariable import DistanceVariable
    >>> var = DistanceVariable(mesh = mesh, value = (-1., -1., -1., -1., 1., 1., 1., 1.))
    >>> var.calcDistanceFunction() #doctest: +LSM
    >>> answer = (-1.75, -1.25, -.75, -0.25, 0.25, 0.75, 1.25, 1.75)
    >>> print var.allclose(answer) #doctest: +LSM
    1

    A 1D test case with very small dimensions.

    >>> dx = 1e-10
    >>> mesh = Grid1D(dx = dx, nx = 8, communicator=serialComm)
    >>> var = DistanceVariable(mesh = mesh, value = (-1., -1., -1., -1., 1., 1., 1., 1.))
    >>> var.calcDistanceFunction() #doctest: +LSM
    >>> answer = numerix.arange(8) * dx - 3.5 * dx
    >>> print var.allclose(answer) #doctest: +LSM
    1

    A 2D test case to test `_calcTrialValue` for a pathological case.

    >>> dx = 1.
    >>> dy = 2.
    >>> from fipy.meshes import Grid2D
    >>> mesh = Grid2D(dx = dx, dy = dy, nx = 2, ny = 3, communicator=serialComm)
    >>> var = DistanceVariable(mesh = mesh, value = (-1., 1., 1., 1., -1., 1.))

    >>> var.calcDistanceFunction() #doctest: +LSM
    >>> vbl = -dx * dy / numerix.sqrt(dx**2 + dy**2) / 2.
    >>> vbr = dx / 2
    >>> vml = dy / 2.
    >>> crossProd = dx * dy
    >>> dsq = dx**2 + dy**2
    >>> top = vbr * dx**2 + vml * dy**2
    >>> sqrt = crossProd**2 *(dsq - (vbr - vml)**2)
    >>> sqrt = numerix.sqrt(max(sqrt, 0))
    >>> vmr = (top + sqrt) / dsq
    >>> answer = (vbl, vbr, vml, vmr, vbl, vbr)
    >>> print var.allclose(answer) #doctest: +LSM
    1

    The `extendVariable` method solves the following equation for a given
    extensionVariable.

    .. math::

       \nabla u \cdot \nabla \phi = 0

    using the fast marching method with an initial condition defined at
    the zero level set.

    >>> from fipy.variables.cellVariable import CellVariable
    >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 2, communicator=serialComm)
    >>> var = DistanceVariable(mesh = mesh, value = (-1., 1., 1., 1.))
    >>> var.calcDistanceFunction() #doctest: +LSM
    >>> extensionVar = CellVariable(mesh = mesh, value = (-1, .5, 2, -1))
    >>> tmp = 1 / numerix.sqrt(2)
    >>> print var.allclose((-tmp / 2, 0.5, 0.5, 0.5 + tmp)) #doctest: +LSM
    1
    >>> var.extendVariable(extensionVar, order=1) #doctest: +LSM
    >>> print extensionVar.allclose((1.25, .5, 2, 1.25)) #doctest: +LSM
    1
    >>> mesh = Grid2D(dx = 1., dy = 1., nx = 3, ny = 3, communicator=serialComm)
    >>> var = DistanceVariable(mesh = mesh, value = (-1., 1., 1.,
    ...                                               1., 1., 1.,
    ...                                               1., 1., 1.))
    >>> var.calcDistanceFunction(order=1) #doctest: +LSM
    >>> extensionVar = CellVariable(mesh = mesh, value = (-1., .5, -1.,
    ...                                                    2., -1., -1.,
    ...                                                   -1., -1., -1.))

    >>> v1 = 0.5 + tmp
    >>> v2 = 1.5
    >>> tmp1 = (v1 + v2) / 2 + numerix.sqrt(2. - (v1 - v2)**2) / 2
    >>> tmp2 = tmp1 + 1 / numerix.sqrt(2)
    >>> print var.allclose((-tmp / 2, 0.5, 1.5, 0.5, 0.5 + tmp,
    ...                      tmp1, 1.5, tmp1, tmp2)) #doctest: +LSM
    1
    >>> answer = (1.25, .5, .5, 2, 1.25, 0.9544, 2, 1.5456, 1.25)
    >>> var.extendVariable(extensionVar, order=1) #doctest: +LSM
    >>> print extensionVar.allclose(answer, rtol = 1e-4) #doctest: +LSM
    1

    Test case for a bug that occurs when initializing the distance
    variable at the interface. Currently it is assumed that adjacent cells
    that are opposite sign neighbors have perpendicular normal vectors. In
    fact the two closest cells could have opposite normals.

    >>> mesh = Grid1D(dx = 1., nx = 3, communicator=serialComm)
    >>> var = DistanceVariable(mesh = mesh, value = (-1., 1., -1.))
    >>> var.calcDistanceFunction() #doctest: +LSM
    >>> print var.allclose((-0.5, 0.5, -0.5)) #doctest: +LSM
    1

    Testing second order. This example failed with Scikit-fmm.

    >>> mesh = Grid2D(dx = 1., dy = 1., nx = 4, ny = 4, communicator=serialComm)
    >>> var = DistanceVariable(mesh = mesh, value = (-1., -1., 1., 1.,
    ...                                               -1., -1., 1., 1.,
    ...                                               1., 1., 1., 1.,
    ...                                               1, 1, 1, 1))
    >>> var.calcDistanceFunction(order=2) #doctest: +LSM
    >>> answer = [-1.30473785, -0.5, 0.5, 1.49923009,
    ...           -0.5, -0.35355339, 0.5, 1.45118446,
    ...            0.5, 0.5, 0.97140452, 1.76215286,
    ...            1.49923009, 1.45118446, 1.76215286, 2.33721352]
    >>> print numerix.allclose(var, answer, rtol=1e-9) #doctest: +LSM
    True

    ** A test for a bug in both LSMLIB and Scikit-fmm **

    The following test gives different result depending on whether
    LSMLIB or Scikit-fmm is used. There is a deeper problem that is
    related to this issue. When a value becomes "known" after
    previously being a "trial" value it updates its neighbors'
    values. In a second order scheme the neighbors one step away also
    need to be updated (if the in between cell is "known" and the far
    cell is a "trial" cell), but are not in either package.  By luck
    (due to trial values having the same value), the values calculated
    in Scikit-fmm for the following example are correct although an
    example that didn't work for Scikit-fmm could also be constructed.

    >>> mesh = Grid2D(dx = 1., dy = 1., nx = 4, ny = 4, communicator=serialComm)
    >>> var = DistanceVariable(mesh = mesh, value = (-1., -1., -1., -1.,
    ...                                               1.,  1., -1., -1.,
    ...                                               1.,  1., -1., -1.,
    ...                                               1.,  1., -1., -1.))
    >>> var.calcDistanceFunction(order=2) #doctest: +LSM
    >>> var.calcDistanceFunction(order=2) #doctest: +LSM
    >>> answer = [-0.5,        -0.58578644, -1.08578644, -1.85136395,
    ...            0.5,         0.29289322, -0.58578644, -1.54389939,
    ...            1.30473785,  0.5,        -0.5,        -1.5,
    ...            1.49547948,  0.5,        -0.5,        -1.5]

    The 3rd and 7th element are different for LSMLIB. This is because
    the 15th element is not "known" when the "trial" value for the 7th
    element is calculated. Scikit-fmm calculates the values in a
    slightly different order so gets a seemingly better answer, but
    this is just chance.

    >>> print numerix.allclose(var, answer, rtol=1e-9) #doctest: +SKFMM
    True

    """
    def __init__(self, mesh, name = '', value = 0., unit = None, hasOld = 0):
        """
        Creates a `distanceVariable` object.

        :Parameters:
          - `mesh`: The mesh that defines the geometry of this variable.
          - `name`: The name of the variable.
	  - `value`: The initial value.
	  - `unit`: the physical units of the variable
          - `hasOld`: Whether the variable maintains an old value.

        """
        CellVariable.__init__(self, mesh, name = name, value = value, unit = unit, hasOld = hasOld)
        self._markStale()

    def _calcValue(self):
        return self._value

    def extendVariable(self, extensionVariable, order=2):
        """

        Calculates the extension of `extensionVariable` from the zero
        level set.

        :Parameters:
          - `extensionVariable`: The variable to extend from the zero
            level set.

        """

        dx, shape = self.getLSMshape()
        extensionValue = numerix.reshape(extensionVariable.value, shape)
        phi = numerix.reshape(self._value, shape)

        if LSM_SOLVER == 'lsmlib':
            from pylsmlib import computeExtensionFields as extension_velocities
        elif LSM_SOLVER == 'skfmm':
            from skfmm import extension_velocities
        else:
            raise Exception, "Neither `lsmlib` nor `skfmm` can be found on the $PATH"

        tmp, extensionValue = extension_velocities(phi, extensionValue, ext_mask=phi < 0., dx=dx, order=order)
        extensionVariable[:] = extensionValue.flatten()

    def getLSMshape(self):
        mesh = self.mesh

        if hasattr(mesh, 'nz'):
            raise Exception, "3D meshes not yet implemented"
        elif hasattr(mesh, 'ny'):
            dx = (mesh.dy, mesh.dx)
            shape = (mesh.ny, mesh.nx)
        elif hasattr(mesh, 'nx'):
            dx = (mesh.dx,)
            shape = mesh.shape
        else:
            raise Exception, "Non grid meshes can not be used for solving the FMM."

        return dx, shape

    def calcDistanceFunction(self, order=2):
        """
        Calculates the `distanceVariable` as a distance function.

        :Parameters:
          - `order`: The order of accuracy for the distance funtion
            calculation, either 1 or 2.

        """

        dx, shape = self.getLSMshape()

        if LSM_SOLVER == 'lsmlib':
            from pylsmlib import distance
        elif LSM_SOLVER == 'skfmm':
            from skfmm import distance
        else:
            raise Exception, "Neither `lsmlib` nor `skfmm` can be found on the $PATH"

        self._value = distance(numerix.reshape(self._value, shape), dx=dx, order=order).flatten()
        self._markFresh()

    @property
    def cellInterfaceAreas(self):
        """
        Returns the length of the interface that crosses the cell

        A simple 1D test:

        >>> from fipy.meshes import Grid1D
        >>> mesh = Grid1D(dx = 1., nx = 4)
        >>> distanceVariable = DistanceVariable(mesh = mesh,
        ...                                     value = (-1.5, -0.5, 0.5, 1.5))
        >>> answer = CellVariable(mesh=mesh, value=(0, 0., 1., 0))
        >>> print numerix.allclose(distanceVariable.cellInterfaceAreas,
        ...                        answer)
        True

        A 2D test case:

        >>> from fipy.meshes import Grid2D
        >>> from fipy.variables.cellVariable import CellVariable
        >>> mesh = Grid2D(dx = 1., dy = 1., nx = 3, ny = 3)
        >>> distanceVariable = DistanceVariable(mesh = mesh,
        ...                                     value = (1.5, 0.5, 1.5,
        ...                                              0.5,-0.5, 0.5,
        ...                                              1.5, 0.5, 1.5))
        >>> answer = CellVariable(mesh=mesh,
        ...                       value=(0, 1, 0, 1, 0, 1, 0, 1, 0))
        >>> print numerix.allclose(distanceVariable.cellInterfaceAreas, answer)
        True

        Another 2D test case:

        >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
        >>> from fipy.variables.cellVariable import CellVariable
        >>> distanceVariable = DistanceVariable(mesh = mesh,
        ...                                     value = (-0.5, 0.5, 0.5, 1.5))
        >>> answer = CellVariable(mesh=mesh,
        ...                       value=(0, numerix.sqrt(2) / 4,  numerix.sqrt(2) / 4, 0))
        >>> print numerix.allclose(distanceVariable.cellInterfaceAreas,
        ...                        answer)
        True

        Test to check that the circumfrence of a circle is, in fact,
        :math:`2\pi r`.

        >>> mesh = Grid2D(dx = 0.05, dy = 0.05, nx = 20, ny = 20)
        >>> r = 0.25
        >>> x, y = mesh.cellCenters
        >>> rad = numerix.sqrt((x - .5)**2 + (y - .5)**2) - r
        >>> distanceVariable = DistanceVariable(mesh = mesh, value = rad)
        >>> print numerix.allclose(distanceVariable.cellInterfaceAreas.sum(), 1.57984690073)
        1


        """
        from fipy.variables.interfaceAreaVariable import _InterfaceAreaVariable
        return _InterfaceAreaVariable(self)

    @property
    def _cellInterfaceNormals(self):
        """

        Returns the interface normals over the cells.

           >>> from fipy.meshes import Grid2D
           >>> from fipy.variables.cellVariable import CellVariable
           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
           >>> distanceVariable = DistanceVariable(mesh = mesh,
           ...                                     value = (-0.5, 0.5, 0.5, 1.5))
           >>> v = 1 / numerix.sqrt(2)
           >>> answer = CellVariable(mesh=mesh,
           ...                       value=(((0, 0, v, 0),
           ...                               (0, 0, 0, 0),
           ...                               (0, 0, 0, 0),
           ...                               (0, v, 0, 0)),
           ...                              ((0, 0, v, 0),
           ...                               (0, 0, 0, 0),
           ...                               (0, 0, 0, 0),
           ...                               (0, v, 0, 0))))
           >>> print numerix.allclose(distanceVariable._cellInterfaceNormals, answer)
           True

        """

        dim = self.mesh.dim

        valueOverFaces = numerix.repeat(self._cellValueOverFaces[numerix.newaxis, ...], dim, axis=0)
        cellFaceIDs = self.mesh.cellFaceIDs
        if cellFaceIDs.shape[-1] > 0:
            interfaceNormals = self._interfaceNormals[...,cellFaceIDs]
        else:
            interfaceNormals = 0

        return MA.where(valueOverFaces < 0, 0, interfaceNormals)

    @property
    def _interfaceNormals(self):
        """

        Returns the normals on the boundary faces only, the other are set to zero.

           >>> from fipy.meshes import Grid2D
           >>> from fipy.variables.faceVariable import FaceVariable
           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
           >>> distanceVariable = DistanceVariable(mesh = mesh,
           ...                                     value = (-0.5, 0.5, 0.5, 1.5))
           >>> v = 1 / numerix.sqrt(2)
           >>> answer = FaceVariable(mesh=mesh,
           ...                       value=((0, 0, v, 0, 0, 0, 0, v, 0, 0, 0, 0),
           ...                              (0, 0, v, 0, 0, 0, 0, v, 0, 0, 0, 0)))
           >>> print numerix.allclose(distanceVariable._interfaceNormals, answer)
           True

        """

        M = self.mesh.dim
        interfaceFlag = numerix.repeat(self._interfaceFlag[numerix.newaxis, ...], M, axis=0)
        return numerix.where(interfaceFlag, self._levelSetNormals, 0)

    @property
    def _interfaceFlag(self):
        """

        Returns 1 for faces on boundary and 0 otherwise.

           >>> from fipy.meshes import Grid2D
           >>> from fipy.variables.faceVariable import FaceVariable
           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
           >>> distanceVariable = DistanceVariable(mesh = mesh,
           ...                                     value = (-0.5, 0.5, 0.5, 1.5))
           >>> answer = FaceVariable(mesh=mesh,
           ...                       value=(0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0))
           >>> print numerix.allclose(distanceVariable._interfaceFlag, answer)
           True

        """
        adjacentCellIDs = self.mesh._adjacentCellIDs
        val0 = numerix.take(numerix.array(self._value), adjacentCellIDs[0])
        val1 = numerix.take(numerix.array(self._value), adjacentCellIDs[1])

        return numerix.where(val1 * val0 < 0, 1, 0)

    @property
    def _cellInterfaceFlag(self):
        """

        Returns 1 for those cells on the interface:

        >>> from fipy.meshes import Grid2D
        >>> from fipy.variables.cellVariable import CellVariable
        >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
        >>> distanceVariable = DistanceVariable(mesh = mesh,
        ...                                     value = (-0.5, 0.5, 0.5, 1.5))
        >>> answer = CellVariable(mesh=mesh, value=(0, 1, 1, 0))
        >>> print numerix.allclose(distanceVariable._cellInterfaceFlag, answer)
        True

        """
        from fipy.variables.interfaceFlagVariable import _InterfaceFlagVariable
        return _InterfaceFlagVariable(self)

    @property
    def _cellValueOverFaces(self):
        """

        Returns the cells values at the faces.

           >>> from fipy.meshes import Grid2D
           >>> from fipy.variables.cellVariable import CellVariable
           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
           >>> distanceVariable = DistanceVariable(mesh = mesh,
           ...                                     value = (-0.5, 0.5, 0.5, 1.5))
           >>> answer = CellVariable(mesh=mesh,
           ...                       value=((-.5, .5, .5, 1.5),
           ...                              (-.5, .5, .5, 1.5),
           ...                              (-.5, .5, .5, 1.5),
           ...                              (-.5, .5, .5, 1.5)))
           >>> print numerix.allclose(distanceVariable._cellValueOverFaces, answer)
           True

        """

        M = self.mesh._maxFacesPerCell
        N = self.mesh.numberOfCells
        return numerix.reshape(numerix.repeat(numerix.array(self._value)[numerix.newaxis, ...], M, axis=0), (M, N))

    @property
    def _levelSetNormals(self):
        """

        Return the face level set normals.

           >>> from fipy.meshes import Grid2D
           >>> from fipy.variables.faceVariable import FaceVariable
           >>> mesh = Grid2D(dx = .5, dy = .5, nx = 2, ny = 2)
           >>> distanceVariable = DistanceVariable(mesh = mesh,
           ...                                     value = (-0.5, 0.5, 0.5, 1.5))
           >>> v = 1 / numerix.sqrt(2)
           >>> answer = FaceVariable(mesh=mesh,
           ...                       value=((0, 0, v, v, 0, 0, 0, v, 0, 0, v, 0),
           ...                              (0, 0, v, v, 0, 0, 0, v, 0, 0, v, 0)))
           >>> print numerix.allclose(distanceVariable._levelSetNormals, answer)
           True
        """

        faceGrad = self.grad.arithmeticFaceValue
        faceGradMag = numerix.array(faceGrad.mag)
        faceGradMag = numerix.where(faceGradMag > 1e-10,
                                    faceGradMag,
                                    1e-10)
        faceGrad = numerix.array(faceGrad)

        ## set faceGrad zero on exteriorFaces
        exteriorFaces = self.mesh.exteriorFaces
        if len(exteriorFaces.value) > 0:
            faceGrad[..., exteriorFaces.value] = 0.

        return faceGrad / faceGradMag

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
