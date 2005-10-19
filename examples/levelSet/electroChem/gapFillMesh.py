#!/usr/bin/env python

"""

The `GapFillMesh` object glues 3 meshes together to form a composite
mesh. The first mesh is a `Grid2D` object that is fine and deals with
the area around the trench or via. The second mesh is a
`GmshImporter2D` object that forms a transition mesh from a fine to a
course region. The third mesh is another `Grid2D` object that forms
the boundary layer. This region consists of very large elements and is
only used for the diffusion in the boundary layer.

"""

__docformat__ = 'restructuredtext'

import Numeric

from fipy.meshes.numMesh.gmshImport import GmshImporter2D
from fipy.meshes.grid2D import Grid2D
from fipy.meshes.numMesh.mesh2D import Mesh2D
import os

class GapFillMesh(Mesh2D):
    """
    
    The following test case tests for diffusion across the domain.

        >>> import Numeric

        >>> from fipy.equations.diffusionEquation import DiffusionEquation
        >>> from fipy.solvers.linearPCGSolver import LinearPCGSolver
        >>> from fipy.boundaryConditions.fixedValue import FixedValue
        >>> from fipy.iterators.iterator import Iterator
        >>> from fipy.variables.cellVariable import CellVariable
        >>> import fipy.tools.dump as dump
        >>> domainHeight = 5.
        
        >>> mesh = GapFillMesh(transitionRegionHeight = 2.,
        ...                    cellSize = 0.1,
        ...                    desiredFineRegionHeight = 1.,
        ...                    desiredDomainHeight = domainHeight,
        ...                    desiredDomainWidth = 1.)
        
        >>> dump.write(mesh, 'tmpmesh')
        >>> mesh = dump.read('tmpmesh')

        >>> mesh.getNumberOfCells() - len(mesh.getCellIDsAboveFineRegion())
        90
        
        >>> var = CellVariable(name = "solution variable",
        ...                    mesh = mesh,
        ...                    value = 0.)
        
        >>> eq =  DiffusionEquation(var,
        ...                         transientCoeff = 0.,
        ...                         diffusionCoeff = 1.,
        ...                         boundaryConditions = (FixedValue(mesh.getBottomFaces(), 0.),
        ...                                               FixedValue(mesh.getTopFaces(), domainHeight)),
        ...                         solver = LinearPCGSolver(tolerance = 1.e-15, 
        ...                                                  steps = 1000))
        
        >>> it = Iterator((eq,))
        >>> it.timestep()

    Evaluate the result:
       
        >>> centers = mesh.getCellCenters()[:,1]
        >>> localErrors = ((centers - Numeric.array(var)) / centers)**2
        >>> globalError = Numeric.sqrt(Numeric.sum(localErrors) / mesh.getNumberOfCells())
        >>> argmax = Numeric.argmax(localErrors)
        >>> print 'maximum local error:',Numeric.sqrt(localErrors[argmax])
        maximum local error: 0.0467527693696
        >>> print 'global error:',globalError
        global error: 0.0171020850885

    """
    
    def __init__(self, cellSize = None, desiredDomainWidth = None, desiredDomainHeight = None, desiredFineRegionHeight = None, transitionRegionHeight = None):

        """

        Arguments:

        `cellSize` - The cell size in the fine grid around the trench.

        `desiredDomainWidth` - The desired domain width.

        `desiredDomainHeight` - The total desired height of the
        domain.

        `desiredFineRegionHeight` - The desired height of the in the
        fine region around the trench.

        `transitionRegionHeight` - The height of the transition
        region.

        """
        self.cellSize = cellSize
        ## Calculate the fine region cell counts.
        nx = int(desiredDomainWidth / self.cellSize)
        ny = int(desiredFineRegionHeight / self.cellSize)

        ## Calculate the actual mesh dimensions
        self.actualFineRegionHeight = ny * self.cellSize
        actualDomainWidth = nx * self.cellSize
        numberOfBoundaryLayerCells = int((desiredDomainHeight - self.actualFineRegionHeight - transitionRegionHeight) / actualDomainWidth)
        self.epsilon = self.cellSize * 1e-10
        self.actualDomainHeight = self.actualFineRegionHeight + transitionRegionHeight + numberOfBoundaryLayerCells * actualDomainWidth
        
        ## Build the fine region mesh.
        self.fineMesh = Grid2D(nx = nx, ny = ny, dx = cellSize, dy = cellSize)

        ## Build the transition mesh and displace.
        transitionMesh = self.buildTransitionMesh(nx, transitionRegionHeight, cellSize) + (0, self.actualFineRegionHeight)

        ## Build the boundary layer mesh.
        boundaryLayerMesh = Grid2D(dx = actualDomainWidth,
                                   dy = actualDomainWidth,
                                   nx = 1,
                                   ny = numberOfBoundaryLayerCells) + (0, self.actualFineRegionHeight + transitionRegionHeight)

        ## Add the meshes together.
        mesh = self.fineMesh._concatenate(transitionMesh, self.epsilon)
        mesh = mesh._concatenate(boundaryLayerMesh, self.epsilon)

        ## Initialize the mesh.
        dict = mesh.__getstate__()
        Mesh2D.__init__(self,dict['vertexCoords'], dict['faceVertexIDs'], dict['cellFaceIDs'])

    def __getstate__(self):
        dict = Mesh2D.__getstate__(self)
        dict['epsilon'] = self.epsilon
        dict['actualDomainHeight'] = self.actualDomainHeight
        dict['actualFineRegionHeight'] = self.actualFineRegionHeight
        dict['cellSize'] = self.cellSize
        dict['fineMesh'] = self.fineMesh
        return dict
            
    def __setstate__(self, dict):
        self.epsilon = dict['epsilon']
        self.actualDomainHeight = dict['actualDomainHeight']
        self.actualFineRegionHeight = dict['actualFineRegionHeight']
        self.cellSize = dict['cellSize']
        self.fineMesh = dict['fineMesh']
        Mesh2D.__init__(self, dict['vertexCoords'], dict['faceVertexIDs'], dict['cellFaceIDs'])
    
    def buildTransitionMesh(self, nx, height, cellSize):

        file = open('mesh.geo', 'w')    

        ## Required to coerce gmsh to use the correct number of
        ## cells in the X direction.
        
        fakeCellSize = nx * cellSize / (nx +0.5)
 
        file.write('cellsize = ' + str(fakeCellSize) + """ ;
        height = """ + str(height) + """ ;
        spacing = """ + str(nx * cellSize) + """ ;
        Point(1) = {0  , 0, 0, cellsize } ;  
        Point(2) = {spacing, 0, 0, cellsize } ;
        Point(3) = {0  , height , 0, spacing } ;
        Point(4) = {spacing, height, 0, spacing } ;
        Line(5) = {1, 2} ;
        Line(6) = {2, 4} ;
        Line(7) = {4, 3} ;
        Line(8) = {3, 1} ;
        Line Loop(9) = {5, 6, 7, 8} ;
        Plane Surface(10) = {9} ; """)
        
        file.close()

        os.system('gmsh mesh.geo -2 -v 0')

        return GmshImporter2D("mesh.msh")

    def getTopFaces(self):
        return self.getFaces(lambda face: face.getCenter()[1] > self.actualDomainHeight - self.epsilon)

    def getBottomFaces(self):
        return self.getFaces(lambda face: face.getCenter()[1] < self.epsilon)

    def getCellIDsAboveFineRegion(self):
        cells = self.getCells(lambda cell: cell.getCenter()[1] > self.actualFineRegionHeight - self.cellSize)
        return [cell.getID() for cell in cells]

    def getFineMesh(self):
        return self.fineMesh
        
class TrenchMesh(GapFillMesh):
    """

    The trench mesh takes the parameters generally used to define a
    trench region and recasts then for the general `GapFillMesh`.

    The following test case tests for diffusion across the domain.

        >>> import Numeric

        >>> from fipy.equations.diffusionEquation import DiffusionEquation
        >>> from fipy.solvers.linearPCGSolver import LinearPCGSolver
        >>> from fipy.boundaryConditions.fixedValue import FixedValue
        >>> from fipy.iterators.iterator import Iterator
        >>> from fipy.variables.cellVariable import CellVariable
        >>> import fipy.tools.dump as dump

        >>> cellSize = 0.05e-6
        >>> trenchDepth = 0.5e-6
        >>> boundaryLayerDepth = 50e-6
        >>> domainHeight = 10 * cellSize + trenchDepth + boundaryLayerDepth
        
        >>> mesh = TrenchMesh(trenchSpacing = 1e-6,
        ...                   cellSize = cellSize,
        ...                   trenchDepth = trenchDepth,
        ...                   boundaryLayerDepth = boundaryLayerDepth,
        ...                   aspectRatio = 1.)

        
        >>> dump.write(mesh, 'tmpmesh')
        >>> mesh = dump.read('tmpmesh')

        >>> mesh.getNumberOfCells() - len(mesh.getElectrolyteCells())        
        150

        >>> var = CellVariable(name = "solution variable",
        ...                    mesh = mesh,
        ...                    value = 0.)
        
        >>> eq =  DiffusionEquation(var,
        ...                         transientCoeff = 0.,
        ...                         diffusionCoeff = 1.,
        ...                         boundaryConditions = (FixedValue(mesh.getBottomFaces(), 0.),
        ...                                               FixedValue(mesh.getTopFaces(), domainHeight)),
        ...                         solver = LinearPCGSolver(tolerance = 1.e-15, 
        ...                                                  steps = 1000))
        
        >>> it = Iterator((eq,))
        >>> it.timestep()

    Evaluate the result:
       
        >>> centers = mesh.getCellCenters()[:,1]
        >>> localErrors = ((centers - Numeric.array(var)) / centers)**2
        >>> globalError = Numeric.sqrt(Numeric.sum(localErrors) / mesh.getNumberOfCells())
        >>> argmax = Numeric.argmax(localErrors)
        >>> print 'maximum local error:',Numeric.sqrt(localErrors[argmax])
        maximum local error: 0.050679531208
        >>> print 'global error:',globalError
        global error: 0.0132375017771

    """

    def __init__(self,
                 trenchDepth = None,
                 trenchSpacing = None,
                 boundaryLayerDepth = None,
                 cellSize = None,
                 aspectRatio = None,
                 angle = 0.,
                 bowWidth = 0.,
                 overBumpRadius = 0.,
                 overBumpWidth = 0.):
        """

        `trenchDepth` - Depth of the trench.

        `trenchSpacing` - The distance between the trenches.

        `boundaryLayerDepth` - The depth of the hydrodynamic boundary
        layer.

        `cellSize` - The cell Size.

        `aspectRatio` - trenchDepth / trenchWidth

        `angle` - The angle for the taper of the trench.

        `bowWidth` - The maximum displacement for any bow in the trench shape.

        `overBumpWidth` - The width of the over bump.

        `overBumpRadius` - The radius of the over bump.
        
        """
        
        self.bowWidth = bowWidth

        self.angle = angle

        self.overBumpRadius = overBumpRadius

        self.overBumpWidth = overBumpWidth

        self.trenchDepth = trenchDepth
        
        self.heightBelowTrench = cellSize * 10.

        self.trenchWidth = self.trenchDepth / aspectRatio

        heightAboveTrench = self.trenchDepth / 1.

        fineRegionHeight = self.heightBelowTrench + self.trenchDepth + heightAboveTrench
        transitionHeight = fineRegionHeight * 3.
        self.domainWidth = trenchSpacing / 2.
        domainHeight = self.heightBelowTrench + self.trenchDepth + boundaryLayerDepth
        
        GapFillMesh.__init__(self,
                             cellSize = cellSize,
                             desiredDomainWidth = self.domainWidth,
                             desiredDomainHeight = domainHeight,
                             desiredFineRegionHeight = fineRegionHeight,
                             transitionRegionHeight = transitionHeight)

    def getElectrolyteCells(self):
        def filter(cell):
            x,y = cell.getCenter()
            Y = (y - (self.heightBelowTrench + self.trenchDepth / 2))

            ## taper
            taper = Numeric.tan(self.angle) * Y

            ## bow
            if abs(self.bowWidth) > 1e-12 and (-self.trenchDepth / 2 < Y < self.trenchDepth / 2):
                param1 = self.trenchDepth**2 / 8 / self.bowWidth
                param2 = self.bowWidth / 2
                bow = -Numeric.sqrt((param1 + param2)**2 - Y**2) + param1 - param2
            else:
                bow = 0

            ## over hang

            
            
            Y = y - (self.heightBelowTrench + self.trenchDepth - self.overBumpRadius)

##            Y = 0
            overBump = 0
            if Y > -self.overBumpRadius:
                if self.overBumpRadius > 1e-12:
                    overBump += self.overBumpWidth / 2 * (1 + Numeric.cos(Y * Numeric.pi / self.overBumpRadius))
                if Y > 0:
                    if Y > self.overBumpRadius:
                        overBump -= self.overBumpRadius
                    else:
                        overBump -= self.overBumpRadius - Numeric.sqrt(self.overBumpRadius**2 - Y**2)

##            print overBump
##            raw_input()
            
            if y > self.trenchDepth + self.heightBelowTrench:
                return 1
            elif y < self.heightBelowTrench:
                return 0
##            elif x < self.domainWidth - self.trenchWidth / 2 - taper + bow + overBump:
            elif x > self.trenchWidth / 2 + taper - bow - overBump:
                return 0
            else:
                return 1

        return self.getCells(filter)

    def __getstate__(self):
        dict = GapFillMesh.__getstate__(self)
        dict['trenchDepth'] = self.trenchDepth
        dict['heightBelowTrench'] = self.heightBelowTrench
        dict['domainWidth'] = self.domainWidth
        dict['trenchWidth'] = self.trenchWidth
        dict['angle'] = self.angle
        dict['bowWidth'] = self.bowWidth
        dict['overBumpWidth'] = self.overBumpWidth
        dict['overBumpRadius'] = self.overBumpRadius
        return dict
            
    def __setstate__(self, dict):
        self.trenchDepth = dict['trenchDepth']
        self.heightBelowTrench = dict['heightBelowTrench']
        self.domainWidth = dict['domainWidth']
        self.trenchWidth = dict['trenchWidth']
        self.angle = dict['angle']
        self.bowWidth = dict['bowWidth']
        self.overBumpWidth = dict['overBumpWidth']
        self.overBumpRadius = dict['overBumpRadius']
        GapFillMesh.__setstate__(self, dict)
    

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 

