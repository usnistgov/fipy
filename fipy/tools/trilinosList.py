import PyTrilinos
from PyTrilinos import Epetra

import numpy

import sys

class _TrilinosArray(object):

    shape = property(lambda self:self._shape,lambda self,s:self.reshape(s),doc="The global shape of the array")
    vLength = property(lambda self:self._vLength,doc="The global length of each vector in the MultiVector")
    epetraMap = property(lambda self:self._map,lambda self,m:self.setMap(m),doc="The Epetra map used to distribute the MultiVector")
    comm = property(lambda self:self._comm,doc="The Epetra communicator used by this _TrilinosArray")
    multiVector = property(lambda self:self._mV,doc="The multiVector wrapped by this _TrilinosArray")
    array = property(lambda self:self._mV.array,doc="The array representing the local _TrilinosArray")
    
    def __init__(self,array=None, vLength = None, shape=None,epetraMap=None,init=True):
        self.takes = {}
        if type(shape) is int:
            shape = (shape,)
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
                    if shape is None:
                        shape = array.shape
                    self._shape = shape
                    size = 1
                    strides = ()
                    for i in shape:
                        strides = (size,)+strides
                        size *= i
                    self._strides = strides
                    self._size = size
                    self._dims = len(shape)
                    inds = numpy.arange(size/shape[-1])
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
                    strides = ()
                    for i in shape:
                        strides = (size,)+strides
                        size *= i
                    self._strides = strides
                    self._size = size
                    self._dims = len(shape)
                    inds = numpy.arange(size/self.shape[-1])
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
                isTril = True
                if shape is None:
                    if epetraMap is None:
                        raise TypeError("_TrilinosArray.__init__() needs a map, shape, or array")
                    else:
                        shape = (epetraMap.NumGlobalElements(),)
                if vLength is None:
                    vLength = shape[-1]
                if epetraMap is None:
                    comm = Epetra.PyComm()
                    epetraMap = Epetra.Map(vLength,0,comm)
                else:
                    comm = epetraMap.Comm()
                    vLength = epetraMap.NumGlobalElements()
                shape = shape[:-1]+(vLength,)
                size = 1
                for i in shape[:-1]:
                    size *= i
                mv = Epetra.MultiVector(epetraMap,size)
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
                        for i in f:
                            for k in i.multiVector:
                                mv[curN] = k
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
                            shape = (epetraMap.NumGlobalElements(),)
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
            size = 1
            strides = ()
            for s in shape:
                strides = (size,)+strides
                size *= s
            self._strides = strides
            self._size = size
            self._vLength = vLength
            self._comm = comm
            self._map = epetraMap
            self._mV = mv
            self._dims = len(shape)
            inds = numpy.arange(size/shape[-1])
            if self._dims != 1:
                self._indices = inds.reshape(shape[:-1])
            else:
                self._indices = inds
            self._reprMV = numpy.arange(vLength)
    
    def take(self,ids,axis=None,mode=None):
        currMap = self._map
        comm = self._comm
        myInds = currMap.MyGlobalElements()
        if axis is None:
            if hasattr(ids,'_mV') and self._dims is 1:
                takes = self.takes
                key = str(ids)
                if takes.has_key(key):
                    throughVec,imp,resMap = takes[key]
                    resVec = Epetra.MultiVector(resMap,1)
                    throughVec.Import(self._mV,imp,Epetra.Insert)
                    resVec += throughVec.flatten()
                    return _TrilinosArray(resVec,shape=(resVec.GlobalLength(),))
                else:
                    resMap = ids._map
                    resVec = Epetra.MultiVector(resMap,1)
                    throughMap = Epetra.Map(-1,ids.array.astype(int)[0],0,comm) #MUST be changed on addition of intvectors
                    throughVec = Epetra.MultiVector(throughMap,self._indices.size)
                    imp = Epetra.Import(throughMap,currMap)
                    throughVec.Import(self._mV,imp,Epetra.Insert)
                    resVec+=throughVec.flatten()
                    takes[key] = (throughVec,imp,resMap)
                    return _TrilinosArray(resVec,shape=(resVec.GlobalLength(),))
            else:
                inds = []
                res = []
                strides = self._strides
                for i in ids:
                    i = int(i)
                    ind = 0
                    for s in strides[:-1]:
                        j = i/s
                        r = i%s
                        ind = ind + j
                        i = r
                    inds.append(ind)
                    res.append(i)
                inds = [(inds[i],res[i],) for i in xrange(len(res)) if res[i] in myInds]
                inds = list(numpy.array(inds).T)
                a = self._mV[inds].array
                m = Epetra.Map(-1,a.size,0,comm)
                mv = Epetra.MultiVector(m,a)
                return _TrilinosArray(mv,shape=mv.GlobalLength())
        else:
            dims = self._dims
            myEls = currMap.MyGlobalElements()
            if axis < 0:
                axis += self._dims
            if axis == self._dims-1:
                if hasattr(ids,'_mV'):
                    takes = self.takes
                    key = str(ids)
                    if takes.has_key(key):
                        throughVec,imp,resMap = takes[key]
                        resVec = Epetra.MultiVector(resMap,1)
                        throughVec.Import(self._mV,imp,Epetra.Insert)
                        resVec += throughVec
                        return _TrilinosArray(resVec,shape=self._shape[:-1]+(resVec.GlobalLength(),))
                    else:
                        resMap = ids._map
                        resVec = Epetra.MultiVector(resMap,self._indices.size)
                        throughMap = Epetra.Map(-1,ids.array.astype(int)[0],0,comm) #MUST be changed on addition of intvectors
                        throughVec = Epetra.MultiVector(throughMap,self._indices.size)
                        imp = Epetra.Import(throughMap,currMap)
                        throughVec.Import(self._mV,imp,Epetra.Insert)
                        resVec+=throughVec
                        takes[key] = (throughVec,imp,resMap)
                        return _TrilinosArray(resVec,shape=self._shape[:-1]+(resVec.GlobalLength(),))
                else:
                    els = [currMap.LID(i) for i in ids if ids in myEls]
                    numEls = len(els)
                    a = self._mV[...,els].array
                    m = Epetra.Map(-1,numEls,0,comm)
                    mv = Epetra.MultiVector(m,a)
                    return _TrilinosArray(mv,shape=self._shape[:-1]+(mv.GlobalLength(),))
            inds = self._indices.take(ids.array.astype(int)[0],axis=axis)
            a = self._mV[inds.flatten()].array
            mv = Epetra.MultiVector(currMap,a)
            return _TrilinosArray(mv,shape = inds.shape+(mv.GlobalLength(),))
    
    def fillWith(self,i):
        self._mV.PutScalar(i)
    
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
    
    def setMap(self,m):
        newMV = Epetra.MultiVector(m,self._indices.size)
        imp = Epetra.Import(m,self._map)
        newMV.Import(self._mV,imp,Epetra.Insert)
        self._vLength = m.NumGlobalElements()
        self._reprMV = numpy.arange(self._vLength)
        self._shape = self._shape[:-1]+(self._vLength,)
        self._map = m
        self._mV = newMV
    
    def reshape(self,shape,*args,**copy):
        if copy == {}:
            copy['copy'] = True
        elif not copy.has_key('copy'):
            raise TypeError, "Invalid keyword arguments.  The only valid keyword is 'copy'."
        shape = (args and (shape,)+args) or shape
        if copy['copy']:
            newTril = _TrilinosArray(init=False)
            newTril._mV = self._mV
            newTril._vLength = self._vLength
            newTril._comm = self._comm
            newTril._map = self._map
            newTril._size = self._size
            newTril._reprMV = self._reprMV
            newTril._shape = shape
            newTril._indices = self._indices.reshape(shape[:-1])
            newTril._dims = len(shape)
            return newTril
        else:
            self._shape = shape
            self._indices = self._indices.reshape(shape[:-1])
            self._dims = len(shape)

    def __getslice__(self,i,j):
        return self.__getitem__((slice(i,j,None),))
    
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
                if type(forMV) is int and y is ():
                    res = 0
                    if myInds.__contains__(forMV):
                        res = mv[0,forMV]
                    res = comm.SumAll(res)
                    return res
                elif type(forMV) is _TrilinosArray and forMV.dtype == bool:
                    els = mv[forMV.array.astype('bool')]
                    numEls = els.size
                    m = Epetra.Map(-1,numEls,0,comm)
                    newMV = Epetra.MultiVector(m,els)
                    return _TrilinosArray(newMV)
                else:
                    inds = indices[y]
                    if inds[0] is 0:
                        inds = ()
                    test = test[forMV]
                    if type(test) != numpy.ndarray:
                        test = [test]
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
                    forMV = last
                    y = y[:-1]
                elif y[-1] == Ellipsis:
                    forMV = Ellipsis
                    y = y[:-1]+(Ellipsis,)
                elif y.__contains__(Ellipsis):
                    forMV = last
                    y = y[:-1]
                else:
                    forMV = slice(0, sys.maxint, None)
                if forMV == Ellipsis or (type(forMV) == slice and forMV.stop - forMV.start >= vLength):
                    inds = indices[y]
                    newMV = Epetra.MultiVector(m,mv[inds.reshape(-1)])
                    return _TrilinosArray(newMV,shape=inds.shape+(vLength,))
                else:
                    inds = indices[y]
                    test = test[forMV]
                    if type(test) is not numpy.ndarray:
                        res = []
                        ind = 0
                        if test in myInds:
                            res = mv[inds.reshape(-1)][:,test]
                            ind = comm.MyPID()
                        ind = comm.SumAll(ind)
                        res = comm.GatherAll(res)[ind]
                        res = res.reshape(inds.shape)
                        m = Epetra.Map(res.shape[-1],0,comm)
                        return _TrilinosArray(Epetra.MultiVector(m,res))
                    els = numpy.array([m.LID(k) for k in test if k in myInds])
                    numEls = els.size
                    m = Epetra.Map(-1,numEls,0,comm)
                    newMV = Epetra.MultiVector(m,mv[inds.reshape(-1)][:,els])
                    return _TrilinosArray(newMV,shape=inds.shape+(m.NumGlobalElements(),))

    def __setslice__(self,i,j,a):
        self.__setitem__((slice(i,j,None),),a)
    
    def __setitem__(self,y,a):
        if isinstance(a,_TrilinosArray):
            if a._dims is 1:
                a = a.array[0]
            else:
                a = a.array
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
            test = [i for i in test if i in myInds]
        else:
            indices = indices[y]
        if dims is 1:
            mv[0,test] = a
        else:
            print repr(mv[indices][:,test])
            mv[indices][:,test] = a
            
    def __mul__(self,other):
        res = self.copy()
        if hasattr(other,'_mV'):
            res._mV *= other._mV
            return res
        res._mV *= other
        return res
    
    def __add__(self,other):
        res = self.copy()
        if hasattr(other,'_mV'):
            res._mV += other._mV
            return res
        res._mV += other
        return res

    def __div__(self,other):
        res = self.copy()
        if hasattr(other,'_mV'):
            res._mV /= other._mV
            return res
        res._mV /= other
        return res

    def __sub__(self,other):
        res = self.copy()
        if hasattr(other,'_mV'):
            res._mV -= other._mV
            return res
        res._mV -= other
        return res

    def __rmul__(self,other):
        res = self.copy()
        if hasattr(other,'_mV'):
            res._mV *= other._mV
            return res
        res._mV *= other
        return res
    
    def __radd__(self,other):
        res = self.copy()
        if hasattr(other,'_mV'):
            res._mV += other._mV
            return res
        res._mV += other
        return res
    
    def __rdiv__(self,other):
        if hasattr(other,'_mV'):
            res = other.copy()
            res._mV /= self._mV
            return res
        res = self.__inv__()
        res._mV *= other
        return res

    def __rsub__(self,other):
        if hasattr(other,'_mV'):
            res = other.copy()
            res._mV -= self._mV
            return res
        res = self.__neg__()
        res._mV += other
        return res

    def __imul__(self,other):
        if hasattr(other,'_mV'):
            self._mV *= other._mV
            return self
        self._mV *= other
        return self
    
    def __iadd__(self,other):
        if hasattr(other,'_mV'):
            self._mV += other._mV
            return self
        self._mV += other
        return self
    
    def __idiv__(self,other):
        if hasattr(other,'_mV'):
            self._mV /= other._mV
            return self
        self._mV /= other
        return self

    def __isub__(self,other):
        if hasattr(other,'_mV'):
            self._mV -= other._mV
            return self
        self._mV -= other
        return self
    
    def __inv__(self):
        res = self.copy()
        res._mV = 1/res._mV
        return res

    def __neg__(self):
        res = self.copy()
        res._mV = res._mV.__neg__()
        return res
    
    def __iter__(self):
        if self._dims is 1:
            return self._mV.array[0].__iter__()
        return self._mV.__iter__()
    
    def __str__(self):
        return self.array.reshape(self.shape[:-1]+(-1,)).__str__()

    def __repr__(self):
        return "_TrilinosArray("+self.multiVector.__str__()+" shape = "+self.shape.__str__()+")"

    def __getattr__(self,attr):
        print repr(attr)
        attribute = getattr(self._mV,attr)
        if callable(attribute):
            return lambda *args,**kwargs: _wrap(attribute,self.shape,*args,**kwargs)
        else:
            return attribute

def _wrap(fn,retShape,*args,**kwargs):
    obj = fn(*args,**kwargs)
    if type(obj) == Epetra.MultiVector:
        return _TrilinosArray(obj,shape=retShape)
    return obj

def arange(stop,start = None,step = None,shape=None,map = None):
    if shape is None:
        if step is None:
            if start is None:
                return _TrilinosArray(numpy.arange(stop),epetraMap=map)
            return _TrilinosArray(numpy.arange(stop,start),epetraMap=map)
        return _TrilinosArray(numpy.arange(start,stop,step),epetraMap=map)
    else:
        if step is None:
            if start is None:
                return _TrilinosArray(numpy.arange(stop).reshape(shape[:-1]+(-1,)),epetraMap=map)
            return _TrilinosArray(numpy.arange(stop,start).reshape(shape[:-1]+(-1,)),epetraMap=map)
        return _TrilinosArray(numpy.arange(start,stop,step).reshape(shape[:-1]+(-1,)),epetraMap=map)
