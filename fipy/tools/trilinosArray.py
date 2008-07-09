import PyTrilinos
from PyTrilinos import Epetra

import numpy

IV = 0
V = 1

class trilArr:

    def __init__(self, shp=None, eMap=None, dType='l', \
                 parallel=True, array=None):

        import operator
        if not operator.xor(shp is None, array is None):
            print "FAIL: Must specify either shape or vector."

        if shp is not None:
            self.shape = trilShape(shp)
        elif array is not None:
            self.shape  = trilShape(numpy.array(array).shape)

        if array is None or str(type(array)).count("Epetra") == 0:
            
            if eMap is not None:
                
                self.comm = eMap.Comm()
                self.eMap = eMap
                self.shape.setMap(self.eMap)

            if eMap is None:

                self.comm = Epetra.PyComm()

                if not parallel:
                    if array is None:
                        self.eMap = None
    
                        if dType=='l':
                            self.vector = Epetra.IntVector(NUMERIX.zeros(shp,dType))
                            self.vtype = IV
    
                        elif dType=='f':
    
                            self.vector = Epetra.Vector(NUMERIX.zeros(shp,dType))
                            self.vtype = V
                    else:
                        tmpArray = numpy.array(array).reshape([-1])
                        if dType=='l':
                            self.vector = Epetra.IntVector(tmpArray)
                            self.vtype = IV
                        elif dType=='f':
                            self.vector = Epetra.Vector(tmpArray)
                            self.vtype = V
                             
                elif parallel:

                    self.eMap = Epetra.Map(self.shape.getSize(),0,self.comm)
                    self.shape.setMap(self.eMap)

            if not hasattr(self, "vector"):
                if array is None:
                    if dType == 'l':
    
                        self.vector = Epetra.IntVector(self.eMap)
                        self.vtype = IV
    
                    if dType == 'f':
    
                        self.vector = Epetra.Vector(self.eMap)
                        self.vtype = V
                else:
                    tmpArray = numpy.array(array).reshape([-1])
                    if dType == 'l':
                        self.vector = Epetra.IntVector(self.eMap,tmpArray)
                        self.vtype = IV
                    if dType == 'f':
                        self.vector = Epetra.Vector(self.eMap,tmpArray)
                        self.vtype = V
                        
            self.dtype = dType

        elif array is not None:
            self.vector = array.copy()
            self.comm = array.Comm()
            self.eMap = array.Map()
            if self.eMap.NumMyElements() != self.eMap.NumGlobalElements():
                self.shape = trilShape(array.size)
                if shape is not None:
                    self.shp.reshape(shp)
            else:
                self.shape = trilShape(self.eMap.NumGlobalElements())
                if shp is not None:
                   self.shape.reshape(shp)
            if isinstance(array, Epetra.IntVector):

                self.vtype = IV
                self.dtype = 'l'
                
            elif isinstance(array, Epetra.Vector):

                self.vtype = V
                self.dtype = 'f'
        self.array = self.vector.array
                

    def fillWith(self, value):
        
        if self.vtype==IV:
            
            self.vector.PutValue(value)
            
        else:
            
            self.vector.PutScalar(value)

    def put(self, ids, values):
        self.insertValues(ids, values)

    def insertValues(self, ids, values):

        # this should operate in accordance with the new shapemap method

        if self.eMap is not None:
            elms = list(self.eMap.MyGlobalElements())
            if type(values) != int:
                values = [v for (i,v) in zip(ids,values) if elms.count(i)>0]
            ids = [self.eMap.LID(i) for i in ids if elms.count(i)>0]
        numpy.put(self.array, ids, values)

    def getValues(self, ids):

        idee = [i for i in ids if self.m.MyGlobalElements().count(i)>0]
        return self.vector[idee]

    def _applyFloatFunction(self, f, optarg=None):

        if optarg is None:
            res = f(self.array)
        else:
            res = f(self.array, optarg.array)    
        v = Epetra.Vector(self.eMap, res)
        return trilArr(array=v)

    def arccos(self):
        return self._applyFloatFunction(numpy.arccos)

    def arccosh(self):
        return self._applyFloatFunction(numpy.arccosh)

    def arcsin(self):
        return self._applyFloatFunction(numpy.arcsin)

    def arcsinh(self):
        return self._applyFloatFunction(numpy.arcsinh)

    def arctan(self):
        return self._applyFloatFunction(numpy.arctan)

    def arctanh(self):
        return self._applyFloatFunction(numpy.arctanh)

    def arctan2(self, other):
        return self._applyFloatFunction(numpy.arctan2, other)

    def cos(self):
        return self._applyFloatFunction(numpy.cos)

    def cosh(self):
        return self._applyFloatFunction(numpy.cosh)

    def tan(self):
        return self._applyFloatFunction(numpy.tan)

    def tanh(self):
        return self._applyFloatFunction(numpy.tanh)

    def log10(self):
        return self._applyFloatFunction(numpy.log10)

    def sin(self):
        return self._applyFloatFunction(numpy.sin)

    def sinh(self):
        return self._applyFloatFunction(numpy.sinh)

    def floor(self):
        return self._applyFloatFunction(numpy.floor)

    def ceil(self):
        return self._applyFloatFunction(numpy.ceil)

    def exp(self):
        return self._applyFloatFunction(numpy.exp)
        
    def log(self):
        return self._applyFloatFunction(numpy.log)
        
    def conjugate(self):
        return self._applyFloatFunction(numpy.conjugate)

    def dot(self, other):
        return self.vector.Dot(other.vector)

    def allequal(self, other):
        return numpy.sum(self.array == other.array) == numpy.size(self.array)

    def allclose(self, other, rtol=1.e-5, atol=1.e-8):
        if self.array.shape != other.array.shape:
            return False
        return sum(1 - (numpy.abs(self.array-other.array) < atol+rtol*numpy.abs(other.array))) == 0

    def globalSum(self):
        return self.comm.SumAll(localSum(self))

    def localSum(self):
        return numpy.sum(self.array)

    def globalTake(self, ids):
        els = list(self.localTake(ids))
        locsize = len(myIDs)
        maxsize = comm.MaxAll(locsize)
        sizes = comm.GatherAll(locsize)       
        while locsize<maxsize:
            els.append(-1)
            locsize=len(els)
        allEls = comm.GatherAll(els)
        allEls = [l for (el,proc) in zip(allEls,range(procs)) for (l,pos) in zip(el,range(sizes[proc]))]
        allEls = numpy.array(allEls).reshape(shp)
        return allEls

    def localTake(self, ids):
        pid = self.comm.MyPID()
        glob = self.eMap.MyGlobalElements()
        indices = numpy.array(indices)
        shp = indices.shape
        indices = indices.reshape(-1)
        myIDs = [m.LID(el) for el in indices if list(glob).count(el)>=1]
        return self[myIDs]

    def reshape(self, shape):
        return self.shape.reshape(shape)

    def getShape(self):
        return self.shape.getShape()

    def getRank(self):
        return self.shape.getRank()

    def __setslice__(self, i, j, y):
        self.__setitem__(slice(i,j,None),y)

    def __getslice__(self, i, j):
        self.__getitem__(slice(i,j,None))
    
    def __setitem__(self, i, y):
        # should operate in accordance with shapemap
        print i,"!",y
        i = self.shape.getLocalIndex(i)
        self.vector.__setitem__(i, y)

    def __getitem__(self, y):
        # should operate in accordance with shapemap
        print y
        y = self.shape.getLocalIndex(y)
        return self.vector.__getitem__(y)

    # needs proper iterator

    def __repr__(self):
        if self.comm.NumProc() == 1:
            return "trilArr("+self._makeArray().__repr__()[6:-1]+")"
        else:
            return "trilArr("+self.vector.array.__repr__()+")"

    def __str__(self):
        if self.comm.NumProc() == 1:
            return self._makeArray().__str__()
        else:
            return self.vector.__str__()

    def _makeArray(self):
	return self.array.reshape(self.shape.getGlobalShape())

    def __or__(self, other):

        return self.array | other.array

class trilShape:

    def __init__(self, shape, eMap=None):
        shape = self._shapeCheck(shape)
        self.globalShape = shape
        self.dimensions = self._dimensions(shape)
        self.actualShape = self._size(shape)
        if eMap is not None:
            self.map = eMap

    def setMap(self, eMap):
        
        if isinstance(eMap,Epetra.Map) or isinstance(eMap,Epetra.BlockMap):
            self.map = eMap
        else:
            print "FAIL: Must be an Epetra Map."

    def getGlobalShape(self):
        return self.globalShape

    def getRank(self):
        return self.dimensions

    def getSize(self):
        return self.actualShape

    def getGlobalIndex(self, index):
        return self._globalTranslateIndices(index)

    def getLocalIndex(self, index):
        ind = self.getGlobalIndex(index)
        return self._globalToLocal(ind)

    def _globalToLocal(self, i):
        if self.map is None:
            return -1
        return self.map.LID(i)

    def _intToSlice(self, i):
        return slice(i,i+1,None)
    
    def _globalTranslateIndices(self, index):

##  ******vaguely borked code, will be fixed ******
##         if type(index)==int:
##             index=[self._intToSlice(index)]
##         elif type(index)==tuple or type(index)==list:
##             while type(index)!=int and len(index)==1:
##                 index=index[0]
##             if type(index)==int:
##                 index=[_intToSlice(index)]
##             else:
                

#        if self._dimensions(index)

        if self._dimensions(index) != self.dimensions:
            return -1
        elif not sum([i<j for (i,j) in zip(index,self.globalShape)]):
            return -2
        
        mult = 1
        lineIndex = 0

        for i in range(len(self.globalShape)+1)[1:]:
            lineIndex += mult*index[-i]
            mult *= self.globalShape[-i]

        return lineIndex

    def _globalTranslateSlice(self, sl):
        
        pass

    def _size(self, shape):
        if type(shape)==tuple or type(shape)==list:
            size = shape[0]
            for i in range(self._dimensions(shape))[1:]:
                size*=shape[i]
        else:
            size = shape
        return size

    def _dimensions(self, shape):
        return len(shape)

    def _shapeCheck(self, shape):
        if type(shape)==int:
            return (shape,)
        if type(shape)==list:
            return tuple(shape)
        if type(shape)!=tuple:
            print "FAIL: Shapes must be ints, lists, or tuples."
        return shape

    def reshape(self, shape):
        shape = self._shapeCheck(shape)
        if self.actualShape != self._size(shape):
            print "FAIL: New shape is differently sized from old shape."
            return -1

        self.globalShape = shape
        self.actualShape = self._size(shape)
        self.dimensions = self._dimensions(shape)

        return 1

    def __str__(self):
        return self.globalShape.__str__()

    def __repr__(self):
        return self.globalShape.__repr__()
    

def isTrilArray(obj):
    return isinstance(obj, trilArr)
