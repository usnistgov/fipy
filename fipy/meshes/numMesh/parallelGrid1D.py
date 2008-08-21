from fipy.tools import numerix
from uniformGrid1D import UniformGrid1D
from PyTrilinos import Epetra
from numerix import MA

class ParallelGrid1D(UniformGrid1D):

    def __init__(self,*args,**kwds):
        self.parallel = True
        UniformGrid1D.__init__(self,*args,**kwds)
        comm = Epetra.PyComm()
        self.cellMap = Epetra.Map(self.numberOfCells,0,comm)
        self.faceMap = Epetra.Map(self.numberOfFaces,0,comm)
        self.exteriorFaces = self.getFacesLeft() | self.getFacesRight()

    
    def _createCells(self):
        self.numberOfFaces = self.nx + 1
        f1 = numerix.arange(self.nx,map=self.cellMap)
        f2 = f1 + 1
        a = numerix.array((f1,f2),cTril = True)
        return a

        
    def _translate(self, vector):
        return ParallelGrid1D(dx = self.dx, nx = self.nx, origin = self.origin + vector)

    def __mul__(self, factor):
        return ParallelGrid1D(dx = self.dx * factor, nx = self.nx, origin = self.origin * factor)

    def getInteriorFaces(self):
        from fipy.variables.faceVariable import FaceVariable
        interiorFaces = FaceVariable(mesh=self, value=False)
        interiorFaces[numerix.arange(self.numberOfFaces-2) + 1] = True
        return interiorFaces
            
    def _getCellIDs(self):
        return numerix.arange(self.numberOfCells,map=self.cellMap)
    
    def _getFaceIDs(self):
        return numerix.arange(self.numberOfFaces,map=self.faceMap)
    
    def _getCellFaceOrientations(self):
        orientations = numerix.ones((2, self.numberOfCells),map=self.cellMap)
        orientations[0] *= -1
        orientations[0,0] = 1
        return orientations

    def _getAdjacentCellIDs(self):
        c1 = numerix.arange(self.numberOfFaces,map=self.faceMap)
        ids = numerix.array((c1 - 1, c1))
        ids[0,0] = ids[1,0]
        ids[1,-1] = ids[0,-1]
        return ids[0], ids[1]

    def _getCellToCellIDs(self):
        c1 = numerix.arange(self.numberOfCells,map=self.cellMap)
        ids = MA.array((c1 - 1, c1 + 1))
        ids[0,0] = MA.masked
        ids[1,-1] = MA.masked
        return ids
        
    def getFaceCellIDs(self):
        c1 = numerix.arange(self.numberOfFaces,map=self.faceMap)
        ids = MA.array((c1-1,c1))
        ids[0,0] = ids[1,0]
        ids[1,0] = MA.masked
        ids[1,-1] = MA.masked
        return ids

    def _getNumberOfLocalFaces(self):
        return self.faceMap.NumMyElements()

    def _getNumberOfLocalCells(self):
        return self.cellMap.NumMyElements()

##     get geometry methods
        
##         from common/mesh
        
    def _getFaceAreas(self):
        return numerix.ones(self.numberOfFaces,'d',map=self.faceMap)

    def _getFaceNormals(self):
        faceNormals = numerix.ones((1, self.numberOfFaces), 'd',map=self.faceMap)
        # The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
        faceNormals[...,0] *= -1
        return faceNormals
    
    def getCellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd',map=self.cellMap) * self.dx

    def getCellCenters(self):
        return ((numerix.arange(self.numberOfCells,map=self.cellMap)[numerix.NewAxis, ...] + 0.5) * self.dx + self.origin) * self.scale['length']

    def _getCellDistances(self):
        distances = numerix.zeros(self.numberOfFaces, 'd',map=self.faceMap)
        distances[1:-1] = self.dx
        distances[0] = self.dx / 2.
        distances[-1] = self.dx / 2.
        return distances

    def _getFaceToCellDistanceRatio(self):
        distances = numerix.ones(self.numberOfFaces, 'd',map=self.faceMap)
        distances *= 0.5
        distances[0] = 1
        distances[-1] = 1
        return distances
    
    def _getFaceTangents1(self):
        return numerix.zeros(self.numberOfFaces, 'd',map=self.faceMap)[numerix.NewAxis, ...]

    def _getFaceTangents2(self):
        return numerix.zeros(self.numberOfFaces, 'd',map=self.faceMap)[numerix.NewAxis, ...]

    def _getCellNormals(self):
        normals = numerix.ones((1, 2, self.numberOfCells), 'd',map=self.cellMap)
        normals[:,0] = -1
        return normals
        
    def _getCellAreas(self):
        return numerix.ones((2, self.numberOfCells), 'd',map=self.cellMap)

    def getFaceCenters(self):
        return numerix.arange(self.numberOfFaces,map=self.faceMap)[numerix.NewAxis, ...] * self.dx + self.origin

    def _getCellVertexIDs(self):
        c1 = numerix.arange(self.numberOfCells,map=self.cellMap)
        return numerix.array((c1 + 1, c1))
