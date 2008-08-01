import PyTrilinos
from PyTrilinos import Epetra

import numpy

import sys

class _TrilinosArray:

    shape = property(lambda self:self._shape,lambda self,s:self.reshape(s,copy=False),"The global shape of the array")
    vLength = property(lambda self:self._vLength,"The global length of each vector in the MultiVector")
    epetraMap = property(lambda self:self._map,lambda self,m:self.setMap(m),"The Epetra map used to distribute the MultiVector")
    comm = property(lambda self:self._comm,"The Epetra communicator used by this _TrilinosArray")
    multiVector = property(lambda self:self._mV,"The multiVector wrapped by this _TrilinosArray")
    array = property(lambda self:self._mV.array,"The array representing the local _TrilinosArray")
    
    def __init__(self,array=None, vLength = None, shape=None,epetraMap=None,init=True):
        #print repr(array),repr(vLength),repr(shape),repr(epetraMap),repr(init)
        if init:
            if type(array) is Epetra.MultiVector:
                if epetraMap is None:
                    self._mV = array
                    self._map = array.Map()
                    self._comm = array.Comm()
                    if vLength is None:
                        vLength = array.GlobalLength()
                        self._vLength = vLength
                    else:
                        self._vLength = vLength
                    if shape is None:
                        self._shape = (vLength,)
                    else:
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
                    return
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
                        shape = nshape+tshape
                        tsize = 1
                        for k in tshape[:-1]:
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
        if numNonNone > dims:
            raise IndexError("invalid index")
        indices = self._indices
        vLength = self._vLength
        mv = self._mV
        m = self._map
        myInds = m.MyGlobalElements()
        comm = self._comm
        test = self._reprMV
        if dims == 1:
            if numNonNone is 0 or y[-1] == Ellipsis or y[-1] == slice(None,None,None):
                shape = indices[y].shape[:-1]+(vLength,) #Have to take off the last index, because this is one-dimensional and so it is automatically 1, and meaningless
                return self.reshape(shape)
            else:
                forMV = y[-1]
                y = y[:-1]
                if type(forMV) is int:
                    res = 0
                    if myInds.__contains__(forMV):
                        res = mv[0,forMV]
                    res = comm.SumAll(res)
                    return res
                elif isinstance(forMV,trilIntArr) and forMV.dtype == bool:
                    els = mv[forMV.array.astype('bool')]
                    numEls = els.size
                    m = Epetra.Map(-1,numEls,0,comm)
                    newMV = Epetra.MultiVector(m,els)
                    return _TrilinosArray(newMV)
                else:
                    inds = indices[y]
                    if inds[0] is 0:
                        inds = ()
                    els = numpy.array([m.LID(k) for k in test[forMV] if k in myInds])
                    numEls = els.size
                    m = Epetra.Map(-1,numEls,0,comm)
                    newMV = Epetra.MultiVector(m,mv[0,els])
                    return _TrilinosArray(newMV,shape=inds.shape[:-1]+(m.NumGlobalElements(),))
        else:
            last = y[-1]
            if numNonNone is 0:
                shape = indices[y].shape +(vLength,)
                return self.reshape(shape)
            else:
                if numNonNone == dims:
                    allInt = True
                    for i in y:
                        if not type(i) is int:
                            allInt = False
                            break
                    if allInt:
                        res = 0
                        if myInds.__contains__(last):
                            res = mv[y]
                        res = comm.SumAll(res)
                        return res
                    forMV = y[-1]
                    y = y[:-1]
                elif y[-1] == Ellipsis:
                    forMV = Ellipsis
                    y = y[:-1]+(Ellipsis,)
                elif y.__contains__(Ellipsis):
                    forMV = y[-1]
                    y = y[:-1]
                else:
                    forMV = slice(0, sys.maxint, None)
                if forMV == Ellipsis or (type(forMV) == slice and forMV.stop - forMV.start >= vLength):
                    inds = indices[y]
                    newMV = Epetra.MultiVector(m,mv[inds.reshape(-1)])
                    return _TrilinosArray(newMV,shape=inds.shape+(vLength,))
                else:
                    inds = indices[y]
                    els = numpy.array([m.LID(k) for k in test[forMV] if k in myInds])
                    numEls = els.size
                    m = Epetra.Map(-1,numEls,0,comm)
                    newMV = Epetra.MultiVector(m,mv[inds][:,els])
                    return _TrilinosArray(newMV,shape=inds.shape+(m.NumGlobalElements(),))
    
    def __setitem__(self,y,a):
        if type(y) is not tuple:
            y = (y,)
        m = self._map
        comm = self._comm
        dims = self._dims
        indices = self._indices
        test = self._reprMV
        myInds = self._map.MyGlobalElements()
        mv = self._mV
        if len(y) == dims or y.__contains__(Ellipsis):
            indices = indices[y[:-1]]
            test = test[y[-1]]
            if type(test) is not numpy.ndarray:
                test = numpy.array([test])
            test = [i for i in test if myInds.__contains__(i)]
        else:
            indices = indices[y]
        print repr(indices),repr(test)
        if dims is 1:
            mv[0,test] = a
        else:
            print repr(mv[indices][:,test])
            mv[indices][:,test] = a
    
    def copy(self):
        newMV = Epetra.MultiVector(self._map,self.array.copy())
        newTril = _TrilinosArray(init=False)
        newTril._mV = newMV
        newTril._shape = self._shape
        newTril._vLength = self._vLength
        newTril._comm = self._comm
        newTril._map = self._map
        newTril._size = self._size
        newTril._dims = self._dims
        newTril._indices = self._indices
        newTril._reprMV = self._reprMV
        return newTril

    def reshape(self,shape,*args,**copy):
        if copy == {}: copy['copy'] = True
        shape = (args and (shape,)+args) or shape
        if copy['copy']:
            newTril = self.copy()
            newTril._shape = shape
            if shape[-1] != self._vLength:
                raise ValueError("The final dimension of a _TrilinosArray must remain unchanged")
            newTril._indices = self._indices.reshape(shape[:-1])
            newTril._dims = len(shape)
            return newTril
        else:
            self._shape = shape[:-1]
            if shape[-1] != self._vLength:
                raise ValueError("The final dimension of a _TrilinosArray must remain unchanged")
            self._indices = self._indices.reshape(shape[:-1])
            self._dims = len(shape)
            
        
    def __str__(self):
        return self.array.reshape(self.shape[:-1]+(-1,)).__str__()

    def __repr__(self):
        return "_TrilinosArray("+self.multiVector.__str__()+" shape = "+self.shape.__str__()+")"

class trilIntArr:
    def __init__(self,array=None,shape=None,map=None,dtype=None):
        pass
