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
from fipy.tools.decorators import getsetDeprecated
from fipy.variables.cellVariable import CellVariable

__all__ = ["DistanceVariable"]

class DistanceVariable(CellVariable):
    r"""
    A `DistanceVariable` object calculates :math:`\phi` so it satisfies,

    .. math::
        
       \abs{\nabla \phi} = 1

    using the fast marching method with an initial condition defined by
    the zero level set.

    Currently the solution is first order, This suffices for initial
    conditions with straight edges (e.g. trenches in
    electrodeposition). The method should work for unstructured 2D grids
    but testing on unstructured grids is untested thus far. This is a 2D
    implementation as it stands. Extending to 3D should be relatively
    simple.

    Here we will define a few test cases. Firstly a 1D test case

    >>> from fipy.meshes import Grid1D
    >>> from fipy.tools import serial
    >>> mesh = Grid1D(dx = .5, nx = 8, communicator=serial)
    >>> from distanceVariable import DistanceVariable
    >>> var = DistanceVariable(mesh = mesh, value = (-1, -1, -1, -1, 1, 1, 1, 1))
    >>> var.calcDistanceFunction()
    >>> answer = (-1.75, -1.25, -.75, -0.25, 0.25, 0.75, 1.25, 1.75)
    >>> print var.allclose(answer)
    1

    A 1D test case with very small dimensions.

    >>> dx = 1e-10
    >>> mesh = Grid1D(dx = dx, nx = 8, communicator=serial)
    >>> var = DistanceVariable(mesh = mesh, value = (-1, -1, -1, -1, 1, 1, 1, 1))
    >>> var.calcDistanceFunction()
    >>> answer = numerix.arange(8) * dx - 3.5 * dx
    >>> print var.allclose(answer)
    1

    A 2D test case to test `_calcTrialValue` for a pathological case.

    >>> dx = 1.
    >>> dy = 2.
    >>> from fipy.meshes import Grid2D
    >>> mesh = Grid2D(dx = dx, dy = dy, nx = 2, ny = 3)
    >>> var = DistanceVariable(mesh = mesh, value = (-1, 1, 1, 1, -1, 1))

    >>> var.calcDistanceFunction()
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
    >>> print var.allclose(answer)
    1

    The `extendVariable` method solves the following equation for a given
    extensionVariable.

    .. math::

       \nabla u \cdot \nabla \phi = 0

    using the fast marching method with an initial condition defined at
    the zero level set. Essentially the equation solves a fake distance
    function to march out the velocity from the interface.

    >>> from fipy.variables.cellVariable import CellVariable
    >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 2)
    >>> var = DistanceVariable(mesh = mesh, value = (-1, 1, 1, 1))
    >>> var.calcDistanceFunction()
    >>> extensionVar = CellVariable(mesh = mesh, value = (-1, .5, 2, -1))
    >>> tmp = 1 / numerix.sqrt(2)
    >>> print var.allclose((-tmp / 2, 0.5, 0.5, 0.5 + tmp))
    1
    >>> var.extendVariable(extensionVar)
    >>> print extensionVar.allclose((1.25, .5, 2, 1.25))
    1
    >>> mesh = Grid2D(dx = 1., dy = 1., nx = 3, ny = 3)
    >>> var = DistanceVariable(mesh = mesh, value = (-1, 1, 1,
    ...                                               1, 1, 1,
    ...                                               1, 1, 1))
    >>> var.calcDistanceFunction()
    >>> extensionVar = CellVariable(mesh = mesh, value = (-1, .5, -1,
    ...                                                    2, -1, -1,
    ...                                                   -1, -1, -1))

    >>> v1 = 0.5 + tmp
    >>> v2 = 1.5
    >>> tmp1 = (v1 + v2) / 2 + numerix.sqrt(2. - (v1 - v2)**2) / 2
    >>> tmp2 = tmp1 + 1 / numerix.sqrt(2)
    >>> print var.allclose((-tmp / 2, 0.5, 1.5, 0.5, 0.5 + tmp, 
    ...                      tmp1, 1.5, tmp1, tmp2))
    1
    >>> answer = (1.25, .5, .5, 2, 1.25, 0.9544, 2, 1.5456, 1.25)
    >>> var.extendVariable(extensionVar)
    >>> print extensionVar.allclose(answer, rtol = 1e-4)
    1

    Test case for a bug that occurs when initializing the distance
    variable at the interface. Currently it is assumed that adjacent cells
    that are opposite sign neighbors have perpendicular normal vectors. In
    fact the two closest cells could have opposite normals.

    >>> mesh = Grid1D(dx = 1., nx = 3)
    >>> var = DistanceVariable(mesh = mesh, value = (-1, 1, -1))
    >>> var.calcDistanceFunction()
    >>> print var.allclose((-0.5, 0.5, -0.5))
    1

    For future reference, the minimum distance for the interface cells can
    be calculated with the following functions. The trial cell values will
    also be calculated with these functions. In essence it is not
    difficult to calculate the level set distance function on an
    unstructured 3D grid. However a lot of testing will be required. The
    minimum distance functions will take the following form.

    .. math::

       X_{\text{min}} = \frac{\left| \vec{s} \times \vec{t} \right|} {\left|
       \vec{s} - \vec{t} \right|}

    and in 3D,

    .. math::
        
       X_{\text{min}} = \frac{1}{3!} \left| \vec{s} \cdot \left( \vec{t} \times
       \vec{u} \right) \right|

    where the vectors :math:`\vec{s}`, :math:`\vec{t}` and :math:`\vec{u}` represent the
    vectors from the cell of interest to the neighboring cell.
    """
    def __init__(self, mesh, name = '', value = 0., unit = None, hasOld = 0, narrowBandWidth = 1e+10):
        """
        Creates a `distanceVariable` object.

        :Parameters:
          - `mesh`: The mesh that defines the geometry of this variable.
          - `name`: The name of the variable.
	  - `value`: The initial value.
	  - `unit`: the physical units of the variable
          - `hasOld`: Whether the variable maintains an old value.
          - `narrowBandWidth`: The width of the region about the zero level set
            within which the distance function is evaluated.

        """
        CellVariable.__init__(self, mesh, name = name, value = value, unit = unit, hasOld = hasOld)
        self._markStale()
        self.narrowBandWidth = narrowBandWidth

        self.cellToCellDistances = MA.filled(self.mesh._cellToCellDistances, 0)
        self.cellNormals = MA.filled(self.mesh._cellNormals, 0)      
        self.cellAreas = MA.filled(self.mesh._cellAreas, 0)
##         self.cellToCellDistances = numerix.array(MA.array(self.mesh._cellToCellDistances).filled(0))
##         self.cellNormals = numerix.array(MA.array(self.mesh._cellNormals).filled(0))       
##         self.cellAreas = numerix.array(MA.array(self.mesh._cellAreas).filled(0))
        self.cellToCellIDs = numerix.array(self.mesh._cellToCellIDsFilled)
        self.adjacentCellIDs = self.mesh._adjacentCellIDs
        self.exteriorFaces = self.mesh.exteriorFaces
        self.cellFaceIDs = self.mesh.cellFaceIDs
        
    def _calcValue(self):
        return self._value
        
    def extendVariable(self, extensionVariable, deleteIslands = False):
        """
        
        Takes a `cellVariable` and extends the variable from the zero
        to the region encapuslated by the `narrowBandWidth`.

        :Parameters:
          - `extensionVariable`: The variable to extend from the zero
            level set.
          - `deleteIslands`: Sets the temporary level set value to
            zero in isolated cells.

        """
        
        self.tmpValue = self._value.copy()
        numericExtensionVariable = numerix.array(extensionVariable)
        self._calcDistanceFunction(numericExtensionVariable, deleteIslands = deleteIslands)
        extensionVariable[:] = numericExtensionVariable
        self._value = self.tmpValue

    def calcDistanceFunction(self, narrowBandWidth = None, deleteIslands = False):
        """
        Calculates the `distanceVariable` as a distance function.

        :Parameters:
          - `narrowBandWidth`: The width of the region about the zero level set
            within which the distance function is evaluated.
          - `deleteIslands`: Sets the temporary level set value to
            zero in isolated cells.

        """
        self._calcDistanceFunction(narrowBandWidth = narrowBandWidth, deleteIslands = deleteIslands)
        self._markFresh()
    
    def _calcDistanceFunction(self, extensionVariable = None, narrowBandWidth = None, deleteIslands = False):

        if narrowBandWidth == None:
            narrowBandWidth = self.narrowBandWidth

        ## calculate interface values

        cellToCellIDs = self.mesh._cellToCellIDs

        if deleteIslands:
            adjVals = numerix.take(self._value, cellToCellIDs)
            adjInterfaceValues = MA.masked_array(adjVals, mask = (adjVals * self._value) > 0)
            masksum = numerix.sum(numerix.logical_not(MA.getmask(adjInterfaceValues)), 0)
            tmp = MA.logical_and(masksum == 4, self._value > 0)
            self._value = MA.where(tmp, -1, self._value)

        adjVals = numerix.take(self._value, cellToCellIDs)
        adjInterfaceValues = MA.masked_array(adjVals, mask = (adjVals * self._value) > 0)
        dAP = self.mesh._cellToCellDistances
        distances = abs(self._value * dAP / (self._value - adjInterfaceValues))
        indices = MA.argsort(distances, 0)
        sign = (self._value > 0) * 2 - 1

        s = distances[indices[0], numerix.arange(indices.shape[1])]

        if self.mesh.dim == 2:

            t = distances[indices[1], numerix.arange(indices.shape[1])]
            u = distances[indices[2], numerix.arange(indices.shape[1])]

            if indices.shape[1] > 0:
                ns = self.cellNormals[..., indices[0], numerix.arange(indices.shape[1])]
                nt = self.cellNormals[..., indices[1], numerix.arange(indices.shape[1])]
            else:
                ns = MA.zeros(self.cellNormals.shape[:-1] + (0,), 'l')
                nt = MA.zeros(self.cellNormals.shape[:-1] + (0,), 'l')

            signedDistance = MA.where(MA.getmask(s),
                                      self._value,
                                      MA.where(MA.getmask(t),
                                               sign * s,
                                               MA.where(abs(numerix.dot(ns,nt)) < 0.9,
                                                        sign * s * t / MA.sqrt(s**2 + t**2),
                                                        MA.where(MA.getmask(u),
                                                                 sign * s,
                                                                 sign * s * u / MA.sqrt(s**2 + u**2)
                                                                 )
                                                        )
                                               )
                                      )
        else:
            signedDistance = MA.where(MA.getmask(s),
                                      self._value,
                                      sign * s)
            

        self._value = signedDistance

        ## calculate interface flag
        masksum = numerix.sum(numerix.logical_not(MA.getmask(distances)), 0)
        interfaceFlag = (masksum > 0).astype('l')

        ## spread the extensionVariable to the whole interface
        flag = True
        if extensionVariable is None:
            extensionVariable = numerix.zeros(self.mesh.numberOfCells, 'd')
            flag = False
            
        ext = numerix.zeros(self.mesh.numberOfCells, 'd')

        positiveInterfaceFlag = numerix.where(self._value > 0, interfaceFlag, 0)
        negativeInterfaceIDs = numerix.nonzero(numerix.where(self._value < 0, interfaceFlag, 0))[0]

        for id in negativeInterfaceIDs:
            tmp, extensionVariable[...,id] = self._calcTrialValue(id, positiveInterfaceFlag, extensionVariable)

        if flag:
            self._value = self.tmpValue.copy()

        ## evaluate the trialIDs
        adjInterfaceFlag = numerix.take(interfaceFlag, cellToCellIDs)
        hasAdjInterface = (numerix.sum(MA.filled(adjInterfaceFlag, 0), 0) > 0).astype('l')

        trialFlag = numerix.logical_and(numerix.logical_not(interfaceFlag), hasAdjInterface).astype('l')

        trialIDs = list(numerix.nonzero(trialFlag)[0])
        evaluatedFlag = interfaceFlag


        for id in trialIDs:
            self._value[...,id], extensionVariable[id] = self._calcTrialValue(id, evaluatedFlag, extensionVariable)

        while len(trialIDs):

            id = trialIDs[numerix.argmin(abs(numerix.take(self._value, trialIDs)))]

            if abs(self._value[...,id]) > narrowBandWidth / 2:
                break

            trialIDs.remove(id)
            evaluatedFlag[...,id] = 1


            for adjID in MA.filled(cellToCellIDs[...,id], -1):
                if adjID != -1:
                    if not evaluatedFlag[...,adjID]:
                        self._value[...,adjID], extensionVariable[...,adjID] = self._calcTrialValue(adjID, evaluatedFlag, extensionVariable)
                        if adjID not in trialIDs:
                            trialIDs.append(adjID)

        self._value = numerix.array(self._value)

    def _calcTrialValue(self, id, evaluatedFlag, extensionVariable):
        adjIDs = self.cellToCellIDs[...,id]
        adjEvaluatedFlag = numerix.take(evaluatedFlag, adjIDs)
        adjValues = numerix.take(self._value, adjIDs)
        adjValues = numerix.where(adjEvaluatedFlag, adjValues, 1e+10)
        try:
            indices = numerix.argsort(abs(adjValues))
        except TypeError:
            # numpy 1.1 raises a TypeError when using argsort function
            indices = abs(adjValues).argsort()
        sign = (self._value[id] > 0) * 2 - 1
        d0 = self.cellToCellDistances[indices[0], id]
        v0 = self._value[..., adjIDs[indices[0]]]
        e0 = extensionVariable[..., adjIDs[indices[0]]]

        N = numerix.sum(adjEvaluatedFlag)

        index0 = indices[0]
        index1 = indices[1]
        index2 = indices[self.mesh.dim]
        
        if N > 1:
            n0 = self.cellNormals[..., index0, id]
            n1 = self.cellNormals[..., index1, id]

            if self.mesh.dim == 2:
                cross = (n0[0] * n1[1] - n0[1] * n1[0])
            else:
                cross = 0.0
                
            if abs(cross) < 0.1:
                if N == 2:
                    N = 1
                elif N == 3:
                    index1 = index2

        if N == 0:
            raise Exception 
        elif N == 1:
            return v0 + sign * d0, e0
        else:
            d1 = self.cellToCellDistances[index1, id]
            n0 = self.cellNormals[..., index0, id]
            n1 = self.cellNormals[..., index1, id]
            v1 = self._value[..., adjIDs[index1]]
            
            crossProd = d0 * d1 * (n0[0] * n1[1] - n0[1] * n1[0])
            dotProd = d0 * d1 * numerix.dot(n0, n1)
            dsq = d0**2 + d1**2 - 2 * dotProd
            
            top = -v0 * (dotProd - d1**2) - v1 * (dotProd - d0**2)
            sqrt = crossProd**2 *(dsq - (v0 - v1)**2)
            sqrt = numerix.sqrt(max(sqrt, 0))



            dis = (top + sign * sqrt) / dsq

            ## extension variable

            e1 = extensionVariable[..., adjIDs[index1]]
            a0 = self.cellAreas[index0, id]
            a1 = self.cellAreas[index1, id]

            if self._value[id] > 0:
                phi = max(dis, 0)
            else:
                phi = min(dis, 0)

            n0grad = a0 * abs(v0 - phi) / d0
            n1grad = a1 * abs(v1 - phi) / d1
            
            return dis, (e0 * n0grad + e1 * n1grad) / (n0grad + n1grad)

    @getsetDeprecated
    def getCellInterfaceAreas(self):
        return self.cellInterfaceAreas

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
        normals = numerix.array(MA.filled(self._cellInterfaceNormals, 0))
        areas = numerix.array(MA.filled(self.mesh._cellAreaProjections, 0))
        return CellVariable(mesh=self.mesh, 
                            value=numerix.sum(abs(numerix.dot(normals, areas)), axis=0))

    @getsetDeprecated
    def _getCellInterfaceNormals(self):
        return self._cellInterfaceNormals
    
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

        N = self.mesh.numberOfCells
        M = self.mesh._maxFacesPerCell
        dim = self.mesh.dim

        valueOverFaces = numerix.repeat(self._cellValueOverFaces[numerix.newaxis, ...], dim, axis=0)
        if self.cellFaceIDs.shape[-1] > 0:
            interfaceNormals = self._interfaceNormals[...,self.cellFaceIDs]
        else:
            interfaceNormals = 0
        from fipy.tools.numerix import MA
        return MA.where(valueOverFaces < 0, 0, interfaceNormals)

    @getsetDeprecated
    def _getInterfaceNormals(self):
        return self._interfaceNormals

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

    @getsetDeprecated
    def _getInterfaceFlag(self):
        return self._interfaceFlag

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
        adjacentCellIDs = self.adjacentCellIDs
        val0 = numerix.take(numerix.array(self._value), adjacentCellIDs[0])
        val1 = numerix.take(numerix.array(self._value), adjacentCellIDs[1])
        
        return numerix.where(val1 * val0 < 0, 1, 0)

    @getsetDeprecated
    def _getCellInterfaceFlag(self):
        return self._cellInterfaceFlag

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
        flag = MA.filled(numerix.take(self._interfaceFlag, self.cellFaceIDs), 0)

        flag = numerix.sum(flag, axis=0)
        
        return numerix.where(numerix.logical_and(self._value > 0, flag > 0), 1, 0)

    @getsetDeprecated
    def _getCellValueOverFaces(self):
        return self._cellValueOverFaces

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

    @getsetDeprecated
    def _getLevelSetNormals(self):
        return self._levelSetNormals

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
        exteriorFaces = self.exteriorFaces.value
        if len(self.exteriorFaces.value) > 0:
            faceGrad[..., self.exteriorFaces.value] = 0.
        
        return faceGrad / faceGradMag 

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test()         
    



            
            
        
                
