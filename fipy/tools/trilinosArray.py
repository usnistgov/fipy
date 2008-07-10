#!/usr/bin/env python

__docformat__ = 'restructuredtext'

import PyTrilinos
from PyTrilinos import Epetra

import numpy

IV = 0
V = 1

class trilArr:
    """
    trilArr is a wrapper for a Trilinos vector
    allows most of the functionality of a numpy array
    works in parallel
    printing multidimensional arrays doesn't work in parallel
    """
    def __init__(self, shape=None, map=None, dType='l', \
                 parallel=True, array=None):
        """
        Creates a trilArr

        :Parameters:
          - `shape`:the shape of the array.  If passed with an array, it overrides the shape of the array.
          - `map`:an Epetra.Map or Epetra.BlockMap describing how to split the job between processors
          - `dType`:the type of data in the array.  However, due to Trilinos limitations, double will be converted to float, and anything besides float or double will become long
          - `parallel`:whether or not this array should be parallelized
          - `array`:the array to be used.  Any iterable type is accepted here (numpy.array, list, tuple).  Default is all zeros
        """
        import operator
        if shape is None and array is None:
            print "FAIL: Must specify either shape or vector."

        if shape is not None:
            self.shp = trilShape(shape)
        elif array is not None:
            self.shp  = trilShape(numpy.array(array).shape)

        if array is None or str(type(array)).count("Epetra") == 0:
            
            if map is not None:
                
                self.comm = map.Comm()
                self.eMap = map
                self.shp.setMap(self.eMap)

            if map is None:

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
                    self.eMap = Epetra.Map(self.shp.getSize(),0,self.comm)
                    self.shp.setMap(self.eMap)

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
                    mine = self.eMap.MyGlobalElements()
                    mini = min(mine)
                    maxi = max(mine)+1
                    if dType == 'l':
                        self.vector = Epetra.IntVector(self.eMap,tmpArray[mini:maxi])
                        self.vtype = IV
                    if dType == 'f':
                        self.vector = Epetra.Vector(self.eMap,tmpArray[mini:maxi])
                        self.vtype = V
                        
            self.dtype = dType

        elif array is not None:
            self.vector = array.copy()
            self.comm = array.Comm()
            self.eMap = array.Map()
            self.shp.setMap(self.eMap)
            if self.eMap.NumMyElements() != self.eMap.NumGlobalElements():
                self.shp = trilShape(array.size)
                if shape is not None:
                    self.shp.reshape(shape)
            else:
                self.shape = trilShape(self.eMap.NumGlobalElements())
                if shape is not None:
                   self.shp.reshape(shape)
            if isinstance(array, Epetra.IntVector):

                self.vtype = IV
                self.dtype = 'l'
                
            elif isinstance(array, Epetra.Vector):

                self.vtype = V
                self.dtype = 'f'
        self.array = self.vector.array
                

    def fillWith(self, value):
        """
        Fills the matrix with a single value

        :Parameters:
          - `value`:what to fill the array with
        
            >>> t = trilArr(shape=(4,))
            >>> t.fillWith(9)
            >>> t.allElems()
            trilArr([9, 9, 9, 9])
        """
        if self.vtype==IV:
            
            self.vector.PutValue(value)
            
        else:
            
            self.vector.PutScalar(value)

    def put(self, ids, values):
        """
        Puts values into the array

        :Parameters:
          - `ids`: Where to put in the values
          - `values`: The values to put in.  If there are less than there are ids, loops through the list multiple times

            >>> t = trilArr(shape=(4,))
            >>> t.put([0],[5])
            >>> t.allElems()
            trilArr([5, 0, 0, 0])
            >>> t.put([1,2,3],[7,8])
            >>> t.allElems()
            trilArr([5, 7, 8, 7])
        """
        self.insertValues(ids, values)

    def insertValues(self, ids, values):
        """
        Puts values into the array

        :Parameters:
          - `ids`: Where to put in the values
          - `values`: The values to put in.  If there are less than there are ids, loops through the list multiple times

            >>> t = trilArr(shape=(4,))
            >>> t.put([0],[5])
            >>> t.allElems()
            trilArr([5, 0, 0, 0])
            >>> t.put([1,2,3],[7,8])
            >>> t.allElems()
            trilArr([5, 7, 8, 7])
        """
        if self.eMap is not None:
            elms = list(self.eMap.MyGlobalElements())
            if type(values) != int:
                values = [v for (i,v) in zip(ids,values) if elms.count(i)>0]
            ids = [self.eMap.LID(i) for i in ids if list(elms).count(i)>0]
        numpy.put(self.array, ids, values)

    def take(self,ids):
        """
        Takes values out of the array
        
        :Parameters:
          - `ids`: What values to take
        
            >>> t = trilArr([1,2,3,4])
            >>> t.take([1,2])
            [2,3]
         """
        self.globalTake(ids)

    def _applyFloatFunction(self, f, optarg=None):
        """
        Applys a fuunction (with at most one additional argument) to this array and returns it.
        
        :Parameters:
          - `f`: the function to apply
          - `optarg`: an additional argument to the function
        """

        if optarg is None:
            res = f(self.array)
        else:
            res = f(self.array, optarg.array)    
        v = Epetra.Vector(self.eMap, res)
        return trilArr(array=v,shape=self.shp.getGlobalShape())

    def arccos(self):
        """
        arccos of this array

            >>> t = trilArr(array=[1, 1, 1, 1])
            >>> t.arccos()
            trilArr([ 0.,  0.,  0.,  0.])
        """
        return self._applyFloatFunction(numpy.arccos)

    def arccosh(self):
        """
        arccosh of this array

            >>> t = trilArr(array=[1, 1, 1, 1])
            >>> t.arccosh()
            trilArr([ 0.,  0.,  0.,  0.])
        """
        return self._applyFloatFunction(numpy.arccosh)

    def arcsin(self):
        """
        arccos of this array

            >>> t = trilArr(array=[1, 1, 1, 1])
            >>> t.arcsin()
            trilArr([ 1.57079633,  1.57079633,  1.57079633,  1.57079633])
        """
        return self._applyFloatFunction(numpy.arcsin)

    def arcsinh(self):
        """
        arcsinh of this array

            >>> t = trilArr(array=[1, 1, 1, 1])
            >>> t.arcsinh()
            trilArr([ 0.88137359,  0.88137359,  0.88137359,  0.88137359])
        """
        return self._applyFloatFunction(numpy.arcsinh)

    def arctan(self):
        """
        arctan of this array

            >>> t = trilArr(array=[1, 1, 1, 1])
            >>> t.arctan()
            trilArr([ 0.78539816,  0.78539816,  0.78539816,  0.78539816])
        """
        return self._applyFloatFunction(numpy.arctan)

    def arctanh(self):
        """
        arctanh of this array

            >>> t = trilArr(array=[1, 1, 1, 1])
            >>> t.arctanh()
            trilArr([ Inf,  Inf,  Inf,  Inf])
        """
        return self._applyFloatFunction(numpy.arctanh)

    def arctan2(self, other):
        """
        arctan of this array/other

        :Parameters:
          - `other`: The array in the denominator

            >>> n = trilArr(array=[0, 0, 0, 0])
            >>> d = trilArr(array=[1, 1, 1, 1])
            >>> n.arctan2(d)
            trilArr([ 0.,  0.,  0.,  0.])
            >>> d.arctan2(n)
            trilArr([ 1.57079633,  1.57079633,  1.57079633,  1.57079633])
        """
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
        els = self.localTake(ids)
        shape = numpy.array(ids).shape
        if els is None:
            els == []
        els = type(els) == numpy.int32 and [els] or list(els)
        locsize = len(els)
        maxsize = self.comm.MaxAll(locsize)
        sizes = self.comm.GatherAll(locsize)
        procs = self.comm.NumProc()
        while locsize<maxsize:
            els.append(-1)
            locsize=len(els)
        allEls = self.comm.GatherAll(els)
        allEls = [l for (el,proc) in zip(allEls,range(procs)) for (l,pos) in zip(el,range(sizes[proc]))]
        allEls = numpy.array(allEls).reshape(shape)
        return allEls

    def localTake(self, ids):
        pid = self.comm.MyPID()
        glob = self.eMap.MyGlobalElements()
        indices = numpy.array(ids)
        indices = indices.reshape(-1)
        myIDs = [el for el in indices if list(glob).count(el)>=1]
        if myIDs == []: return []
        print "My IDs",self[myIDs]
        return self[myIDs]

    def reshape(self, shape, copy=False):
        ## reshape checks need to be done
        ## before a copy is made
        if self.shp._shapeCheck(shape) is None:
            return
        if copy:
            newArr = self.__copy__()
            shp = newArr.shp
        else:
            shp = self.shp

        shp.reshape(shape)

    def getShape(self):
        return self.shp.getShape()

    def getRank(self):
        return self.shp.getRank()

    def allElems(self):
        """
        Returns the full array
        
            >>> t = trilArr(shape=(4,))
            >>> t.allElems()
            trilArr([0, 0, 0, 0])
        """
        comm = self.vector.Comm()
        pid = comm.MyPID()
        procs = comm.NumProc()
        m = self.vector.Map()
        sz = m.NumGlobalElements()
        locsize = self.vector.MyLength()
        maxsize = comm.MaxAll(locsize)
        els = list(self.vector)
        while locsize<maxsize:
            els.append(-1)
            locsize+=1
        allEls = comm.GatherAll(els)
        allEls = numpy.array(allEls).reshape(-1)
        if sz%procs:
            allEls = [i for (i,j) in zip(allEls,range(1,len(allEls)+1)) \
                      if j<=maxsize*(sz%procs) or j%maxsize]
        return trilArr(array=numpy.array(allEls),shape=self.shp.getGlobalShape(),parallel=False)

    def __setslice__(self, i, j, y):
        self.__setitem__(slice(i,j,None),y)

    def __getslice__(self, i, j):
        self.__getitem__(slice(i,j,None))
    
    def __setitem__(self, i, y):
        # should operate in accordance with shapemap
        i = self.shp.getLocalIndex(i)
        self.vector.__setitem__(i, y)

    def __getitem__(self, y):
        # should operate in accordance with shapemap
        y = self.shp.getLocalIndex(y)
        return self.vector.__getitem__(y)

    def __copy__(self):
        pass

    # needs proper iterator

    def __repr__(self):
        if self.comm.NumProc() == 1:
            return "trilArr("+self._makeArray().__repr__()[6:-1]+")"
        else:
            return "trilArr("+self.vector.array.__repr__()[6:-1]+")"

    def __str__(self):
        if self.comm.NumProc() == 1:
            return self._makeArray().__str__()
        else:
            return self.vector.__str__()

    def _makeArray(self):
	return self.array.reshape(self.shp.getGlobalShape())

    def __or__(self, other):

        return self.array | other.array

class trilShape:

    def __init__(self, shape, eMap=None):
        if str(type(shape)).count("int") != 0: shape = (shape,)
        self.globalShape = shape
        self.dimensions = self._dimensions(shape)
        self.actualShape = self._size(shape)
        shape = self._shapeCheck(shape)
        if eMap is not None:
            self.map = eMap
        mult = 1
        tmp = []
        for i in range(len(self.globalShape)+1)[1:]:
            tmp.append(mult)
            mult *= self.globalShape[-i]
        self.steps = tuple(tmp)

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
    
    def getSteps(self):
        return self.steps
    
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
        if type(i)==slice:
            return i
        else:
            return slice(i,i+1,None)

    def _fillToDim(self, i):
        i = list[i]
        while len(i)<self.dimensions:
            i.append(slice(None,None,None))
        return tuple(i)
    
    def _globalTranslateIndices(self, index):

        if type(index)==int:
            index=[self._intToSlice(index)]
        elif type(index)==tuple or type(index)==list:
            if type(index[0])!=int and type(index[0])!=slice:
                while type(index)!=int and len(index)==1:
                    index=index[0]
                if type(index)==int:
                    index=[self._intToSlice(index)]
                elif len(index)<=self.dimensions:
                    index = [[i[el] for i in index] for el in range(len(i))]
            else:
                index = [tuple(index)]
        index = [self._fillToDim(i) for i in index]

        index = [el for el in self._globalTranslateSlices(i) for i in index]

        indices = []
        for ind in index:
            if self._dimensions(ind)!=self.dimensions:
                return -1
            if not sum([i<j for (i,j) in zip(ind,self.globalShape)]):
                return -2
        
            lineIndex = 0
            for mult in self.steps:
                lineIndex += mult*index[-i]

            indices.append(lineIndex)

        return indices

    def _globalTranslateSlices(self, sls):
        for (el,i) in zip(sls,range(len(sls))):
            if type(el)==int:
                sls[i]=self._intToSlice(el)
        dims = [range(i) for i in self.globalShape]
        res = [tuple(dim[sl]) for (sl,dim) in zip(sls,dims)]
        ## bork bork bork ##

    def _size(self, shape):
        if type(shape)==tuple or type(shape)==list:
            size = shape[0]
            for i in range(self._dimensions(shape))[1:]:
                size*=shape[i]
        else:
            size = shape
        return size

    def _dimensions(self, shape):
        if str(type(shape)).count("int") == 1:
            return 1
        return len(shape)

    def _shapeCheck(self, shape):
        if type(shape)==int:
            shape = (shape,)
        if type(shape)==list:
            shape = tuple(shape)
        if type(shape)!=tuple:
            print "FAIL: Shapes must be ints, lists, or tuples."
            return None
        if self.actualShape != self._size(shape):
            print "FAIL: New shape is differently sized from old shape."
            return None
        return shape

    def reshape(self, shape):
        shape = self._shapeCheck(shape)
        if shape is None:
            return -1

        self.globalShape = shape
        self.actualShape = self._size(shape)
        self.dimensions = self._dimensions(shape)

        mult = 1
        tmp = []
        for i in range(len(self.globalShape)+1)[1:]:
            tmp.append(mult)
            mult *= self.globalShape[-i]
        self.steps = tuple(tmp)

        return 1

    def __str__(self):
        return self.globalShape.__str__()

    def __repr__(self):
        return "trilShape("+self.globalShape.__repr__()+")"
    

def isTrilArray(obj):
    return isinstance(obj, trilArr)

if __name__ == '__main__':
    import doctest
    doctest.testmod()