import PyTrilinos
from PyTrilinos import Epetra

import numpy

class _TrilinosArray:

    shape = property(lambda self:self._shape,"The global shape of the array")
    vLength = property(lambda self:self._vLength,"The global length of each vector in the MultiVector")
    epetraMap = property(lambda self:self._map,lambda self,m:self.setMap(m),"The Epetra map used to distribute the MultiVector")
    comm = property(lambda self:self._comm,"The Epetra communicator used by this _TrilinosArray")
    multiVector = property(lambda self:self._mV,"The multiVector wrapped by this _TrilinosArray")
    array = property(lambda self:self.__array__(),"The array representing the global _TrilinosArray")
    
    def __init__(self,array=None, vLength = None, shape=None,epetraMap=None):
        if type(array) is Epetra.MultiVector:
            if epetraMap is None:
                self._mV = array
                self._map = array.Map()
                self._comm = array.Comm()
                self._vLength = vLength
                self._shape = shape
                size = 1
                for i in shape[:-1]:
                    size *= i
                self._size = size*shape[-1]
                self._dims = len(shape)
                inds = numpy.arange(size)
                if self._dims != 1:
                    self._indices = inds.reshape(shape[:-1])
                else:
                    self._indices = inds
                self._reprMV = numpy.arange(vLength)
            else:
                self._map = epetraMap
                self._comm = array.Comm()
                self._vLength = array.GlobalLength()
                self._shape = shape
                size = 1
                for i in shape[:-1]:
                    size *= i
                self._size = size*shape[-1]
                self._dims = len(shape)
                inds = numpy.arange(size)
                if self._dims != 1:
                    self._indices = inds.reshape(shape[:-1])
                else:
                    self._indices = inds
                self._reprMV = numpy.arange(vLength)
                oldMap = array.Map()
                DistToPers = Epetra.Import(epetraMap,oldMap)
                mv = Epetra.MultiVector(epetraMap,array.NumVectors())
                mv.Import(array, DistToPers, Epetra.Insert)
                self._mV = mv
            return
        if array is None:
            isTril = False
            if shape is None:
                if epetraMap is None:
                    raise TypeError("_TrilinosArray.__init__() needs a map, shape, or array")
                else:
                    shape = (epetraMap.NumGlobalElements())
            if vLength is None:
                vLength = shape[-1]
            if epetraMap is None:
                comm = Epetra.PyComm()
                epetraMap = Epetra.Map(vLength,0,comm)
            else:
                comm = epetraMap.Comm()
                vLength = epetraMap.NumGlobalElements()
                shape = shape[:-1]+(vLength,)
            mv = Epetra.MultiVector(epetraMap)
        else:
            isTril = isinstance(array,_TrilinosArray)
            if not isTril:
                narray = numpy.array(array,dtype='object')
                t1 = narray.take([0])[0]
                isTril = isinstance(t1,_TrilinosArray)
                if isTril:
                    depth = len(narray.shape)
                    nshape = narray.shape
                    tshape = t1.shape
                    tsize = tshape[0]
                    for k in tshape[1:-1]:
                        tsize*=k
                    t1 = t1.multiVector
                    comm = t1.Comm()
                    vLength = t1.GlobalLength()
                    if epetraMap is None:
                        epetraMap = Epetra.Map(vLength,0,comm)
                    mv = Epetra.MultiVector(epetraMap,narray.size*tsize)
                    f = narray.flat
                    curN = 0
                    for i in range(len(f)):
                        v = f[i].multiVector
                        for k in range(len(v)):
                            mv[curN,:] = v[k,:]
                            curN+=1
            else:
                if shape is None or shape == array.shape or shape == (array.shape,):
                    shape = array.shape
                    array = array.multiVector
                    vLength = array.GlobalLength()
                    comm = array.Comm()
                    oldMap = array.Map()
                    if epetraMap is None:
                        epetraMap = oldMap
                        mv = Epetra.MultiVector(Epetra.Copy,array)
                    else:
                        DistToPers = Epetra.Import(epetraMap,oldMap)
                        mv = Epetra.MultiVector(epetraMap,array.NumVectors())
                        mv.Import(array, DistToPers, Epetra.Insert)
                else:
                    array = array.multiVector
                    PersonalMap = Epetra.Map(-1, range(0, array.GlobalLength()), 0, v.Comm())
                    DistToPers = Epetra.Import(PersonalMap, array.Map())
                    PersonalV = Epetra.Vector(PersonalMap)
                    PersonalV.Import(array, DistToPers, Epetra.Insert)
                    vLength = shape[-1]
                    PersonalV = PersonalV.reshape(-1,vLength)
                    comm = array.comm
                    if epetraMap is None:
                        epetraMap = Epetra.Map(vLength,0,comm)
                    else:
                        vLength = epetraMap.NumGlobalElements()
                        shape = shape[:-1]+(vLength,)
                    if type(shape) != tuple:
                        shape = (shape,)
                    PersonalV = PersonalV[...,epetraMap.MyGlobalElements()]
                    mv = Epetra.MultiVector(epetraMap,narray)
        if not isTril:
            del t1
            if shape is None:
                if array is not None:
                    shape = narray.shape
                else:
                    if epetraMap is None:
                        raise TypeError("_TrilinosArray.__init__() needs a map, shape, or array")
                    else:
                        shape = (epetraMap.NumGlobalElements())
            if vLength is None:
                vLength = shape[-1]
            narray = narray.reshape(-1,vLength).astype('float')
            if epetraMap is None:
                comm = Epetra.PyComm()
                epetraMap = Epetra.Map(vLength,0,comm)
            else:
                comm = epetraMap.Comm()
                vLength = epetraMap.NumGlobalElements()
                shape = shape[:-1]+(vLength,)
            narray = narray[...,epetraMap.MyGlobalElements()]
            mv = Epetra.MultiVector(epetraMap,narray)
        self._shape = shape
        self._vLength = vLength
        self._comm = comm
        self._map = epetraMap
        self._mV = mv
        size = 1
        for i in shape[:-1]:
            size *= i
        self._size = size*shape[-1]
        self._dims = len(shape)
        inds = numpy.arange(size)
        if self._dims != 1:
            self._indices = inds.reshape(shape[:-1])
        else:
            self._indices = inds
        self._reprMV = numpy.arange(vLength)
    
    def __getitem__(self,y):
        if type(y) is not tuple:
            y = (y,)
        shape = ()
        numS = len(y)
        numNone = list(y).count(None)
        numNonNone = numS-numNone
        dims = self._dims
        indices = self._indices
        vLength = self._vLength
        mv = self._mV
        m = self._map
        myInds = m.MyGlobalElements()
        comm = self._comm
        test = self._reprMV
        if dims == 1:
            if numNonNone is 0:
                shape = indices[y].shape+(vLength,)
                return self.copy().reshape(shape)
            elif numNone is 0:
                y = y[0]
                if type(y) is int:
                    res = 0
                    if myInds.__contains__(y):
                        res = mv[0,y]
                    res = comm.SumAll(res)
                    return res
                els = test[y]
                allMap = Epetra.Map(-1,els.size,0,comm)
                toDist = Epetra.Import(allMap,m)
                newMV = Epetra.MultiVector(allMap,1)
                newMV.Import(mv,toDist,Epetra.Insert)
                return _TrilinosArray(newMV,vLength = els.size,shape=(els.size,))
        if numNonNone == dims:
            formv = y[-1]
            if formv is None:
                raise ValueError("Once a _TrilinosArray is initialized, you cannot change the vector length")
            y = y[:-1]
        
        
        
    def __str__(self):
        return self.multiVector.__str__()

    def __repr__(self):
        return "_TrilinosArray("+self.multiVector.__str__()+" shape = "+self.shape.__str__()+")"

class trilIntArr:
    def __init__(self,array=None,shape=None,map=None,dtype=None):
        pass
