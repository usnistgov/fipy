#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "setup.py"
 #
 #  Author: Jonathan Guyer   <guyer@nist.gov>
 #  Author: Daniel Wheeler   <daniel.wheeler@nist.gov>
 #  Author: James Warren     <jwarren@nist.gov>
 #  Author: Andrew Acquaviva <andrewa@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  setup.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 # ###################################################################
 ##

"""

The `gapFillMesh` function glues 3 meshes together to form a composite
mesh. The first mesh is a `Grid2D` object that is fine and deals with
the area around the trench or via. The second mesh is a `Gmsh2D`
object that forms a transition mesh from a fine to a course
region. The third mesh is another `Grid2D` object that forms the
boundary layer. This region consists of very large elements and is
only used for the diffusion in the boundary layer.

"""

__docformat__ = 'restructuredtext'

from fipy.meshes import Gmsh2D
from fipy.meshes import Grid2D
from fipy.tools import numerix
from fipy.tools import serial

def gapFillMesh(cellSize=None,
                desiredDomainWidth=None,
                desiredDomainHeight=None,
                desiredFineRegionHeight=None,
                transitionRegionHeight=
                None):
    """

    Arguments:

    `cellSize` - The cell size in the fine grid around the trench.

    `desiredDomainWidth` - The desired domain width.

    `desiredDomainHeight` - The total desired height of the domain.

    `desiredFineRegionHeight` - The desired height of the in the fine
    region around the trench.

    `transitionRegionHeight` - The height of the transition region.

    The following test case tests for diffusion across the domain.

    >>> domainHeight = 5.        
    >>> mesh = gapFillMesh(transitionRegionHeight = 2.,
    ...                    cellSize = 0.1,
    ...                    desiredFineRegionHeight = 1.,
    ...                    desiredDomainHeight = domainHeight,
    ...                    desiredDomainWidth = 1.) # doctest: +GMSH

    >>> import fipy.tools.dump as dump
    >>> (f, filename) = dump.write((mesh, mesh.cellIDsAboveFineRegion)) # doctest: +GMSH
    >>> mesh, cellIDsAboveFineRegion = dump.read(filename, f) # doctest: +GMSH
    >>> mesh.cellIDsAboveFineRegion = cellIDsAboveFineRegion
    >>> mesh.numberOfCells - len(mesh.cellIDsAboveFineRegion) # doctest: +GMSH
    90

    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(mesh = mesh) # doctest: +GMSH

    >>> from fipy.terms.diffusionTerm import DiffusionTerm
    >>> eq = DiffusionTerm()

    >>> var.constrain(0., mesh.facesBottom) # doctest: +GMSH
    >>> var.constrain(domainHeight, mesh.facesTop) # doctest: +GMSH

    >>> eq.solve(var) # doctest: +GMSH

    Evaluate the result:

    >>> centers = mesh.cellCenters[1].copy() # doctest: +GMSH 

    .. note:: the copy makes the array contiguous for inlining

    >>> localErrors = (centers - var)**2 / centers**2 # doctest: +GMSH 
    >>> globalError = numerix.sqrt(numerix.sum(localErrors) / mesh.numberOfCells) # doctest: +GMSH 
    >>> argmax = numerix.argmax(localErrors) # doctest: +GMSH 
    >>> print numerix.sqrt(localErrors[argmax]) < 0.1 # doctest: +GMSH 
    1
    >>> print globalError < 0.05 # doctest: +GMSH 
    1

    """

    # Calculate the fine region cell counts.
    nx = int(desiredDomainWidth / cellSize)
    ny = int(desiredFineRegionHeight / cellSize)

    # Calculate the actual mesh dimensions
    actualFineRegionHeight = ny * cellSize
    actualDomainWidth = nx * cellSize
    numberOfBoundaryLayerCells = int((desiredDomainHeight - actualFineRegionHeight - transitionRegionHeight) / actualDomainWidth)

    # Build the fine region mesh.
    fineMesh = Grid2D(nx = nx, ny = ny, dx = cellSize, dy = cellSize, communicator=serial)

    # Build the transition mesh and displace.
    transitionMesh = buildTransitionMesh(nx, transitionRegionHeight, cellSize) + ((0,), (actualFineRegionHeight,))

    # Build the boundary layer mesh.

    boundaryLayerMesh = Grid2D(dx = actualDomainWidth,
                               dy = actualDomainWidth,
                               nx = 1,
                               ny = numberOfBoundaryLayerCells,
                               communicator=serial) + ((0,), (actualFineRegionHeight + transitionRegionHeight,),)

    # Add the meshes together.
    mesh = fineMesh + transitionMesh + boundaryLayerMesh

    mesh.cellIDsAboveFineRegion = numerix.nonzero(mesh.cellCenters[1] > actualFineRegionHeight - cellSize)[0]
    mesh.fineMesh = fineMesh
    return mesh

def buildTransitionMesh(nx, height, cellSize):

    # Required to coerce gmsh to use the correct number of
    # cells in the X direction.

    fakeCellSize = nx * cellSize / (nx +0.5)

    # including extra point with height extraPointHeight due to new gmsh 2.0 issue, see thread
    # http://www.geuz.org/pipermail/gmsh/2007/002465.html

    return Gmsh2D('cellsize = ' + str(fakeCellSize) + """ ;
    height = """ + str(height) + """ ;
    spacing = """ + str(nx * cellSize) + """ ;
    extraPointHeight = """ + str(100. * height) + """ ;
    Point(1) = {0  , 0, 0, cellsize } ;  
    Point(2) = {spacing, 0, 0, cellsize } ;
    Point(3) = {0  , height , 0, spacing } ;
    Point(4) = {spacing, height, 0, spacing } ;
    Point(5) = {0, extraPointHeight, 0, cellsize} ;
    Line(5) = {1, 2} ;
    Line(6) = {2, 4} ;
    Line(7) = {4, 3} ;
    Line(8) = {3, 1} ;
    Line Loop(9) = {5, 6, 7, 8} ;
    Plane Surface(10) = {9} ; """, communicator=serial)

def trenchMesh(trenchDepth=None,
               trenchSpacing=None,
               boundaryLayerDepth=None,
               cellSize=None,
               aspectRatio=None,
               angle=0.)
    """

    `trenchDepth` - Depth of the trench.

    `trenchSpacing` - The distance between the trenches.

    `boundaryLayerDepth` - The depth of the hydrodynamic boundary
    layer.

    `cellSize` - The cell Size.

    `aspectRatio` - trenchDepth / trenchWidth

    `angle` - The angle for the taper of the trench.

    The trench mesh takes the parameters generally used to define a
    trench region and recasts then for the general `gapFillMesh`.

    The following test case tests for diffusion across the domain.

    >>> cellSize = 0.05e-6
    >>> trenchDepth = 0.5e-6
    >>> boundaryLayerDepth = 50e-6
    >>> domainHeight = 10 * cellSize + trenchDepth + boundaryLayerDepth

    >>> mesh = trenchMesh(trenchSpacing = 1e-6,
    ...                   cellSize = cellSize,
    ...                   trenchDepth = trenchDepth,
    ...                   boundaryLayerDepth = boundaryLayerDepth,
    ...                   aspectRatio = 1.) # doctest: +GMSH

    >>> import fipy.tools.dump as dump
    >>> (f, filename) = dump.write((mesh, mesh.electrolyteMask)) # doctest: +GMSH
    >>> mesh, electrolyteMask = dump.read(filename, f) # doctest: +GMSH
    >>> mesh.electrolyteMask = electrolyteMask
    >>> mesh.numberOfCells - len(numerix.nonzero(mesh.electrolyteMask)[0]) # doctest: +GMSH
    150

    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(mesh = mesh, value = 0.) # doctest: +GMSH

    >>> from fipy.terms.diffusionTerm import DiffusionTerm
    >>> eq = DiffusionTerm() # doctest: +GMSH

    >>> var.constrain(0., mesh.facesBottom) # doctest: +GMSH
    >>> var.constrain(domainHeight, mesh.facesTop) # doctest: +GMSH

    >>> eq.solve(var) # doctest: +GMSH

    Evaluate the result:
       
    >>> centers = mesh.cellCenters[1].copy() # doctest: +GMSH

    .. note:: the copy makes the array contiguous for inlining

    >>> localErrors = (centers - var)**2 / centers**2 # doctest: +GMSH
    >>> globalError = numerix.sqrt(numerix.sum(localErrors) / mesh.numberOfCells) # doctest: +GMSH
    >>> argmax = numerix.argmax(localErrors) # doctest: +GMSH
    >>> print numerix.sqrt(localErrors[argmax]) < 0.051 # doctest: +GMSH
    1
    >>> print globalError < 0.02 # doctest: +GMSH
    1

    """

    heightBelowTrench = cellSize * 10.

    heightAboveTrench = trenchDepth / 1.

    fineRegionHeight = heightBelowTrench + trenchDepth + heightAboveTrench
    transitionHeight = fineRegionHeight * 3.
    domainWidth = trenchSpacing / 2.
    domainHeight = heightBelowTrench + trenchDepth + boundaryLayerDepth

    mesh = gapFillMesh(cellSize=cellSize,
                       desiredDomainWidth=domainWidth,
                       desiredDomainHeight=domainHeight,
                       desiredFineRegionHeight=fineRegionHeight,
                       transitionRegionHeight=transitionHeight)

    trenchWidth = trenchDepth / aspectRatio
    
    x, y = mesh.cellCenters
    Y = (y - (heightBelowTrench + trenchDepth / 2))
    taper = numerix.tan(angle) * Y
    mesh.electrolyteMask = numerix.where(y > trenchDepth + heightBelowTrench,
                                         1,
                                         numerix.where(y < heightBelowTrench,
                                                       0,
                                                       numerix.where(x > trenchWidth / 2 + taper,
                                                                     0,
                                                                     1)))

    return mesh

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test() 

