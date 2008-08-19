import PyTrilinos
from PyTrilinos import Epetra

import numpy

import sys

num = 0
count = False
debug = False

class _TrilinosArray(object):
    
    def __init__(self,array=None, vLength = None, shape=None,map=None,dtype=None,init=True):
        self.takes = {}
        if debug:
            print "IN __init__\n"
        if init:
            if fixType(type(shape)) == 'int':
                shape = (shape,)
            elif type(shape) is list:
                shape = tuple(shape)
            if count:
                global num
                num += 1
                print num,repr(array),vLength,shape,repr(map),dtype
            if array is None and shape is None and map is None:
                raise TypeError("_TrilinosArray.__init__() needs a map, shape, or array")
            if hasattr(array,'getValue'):
                array = array.getValue()
            if dtype is None:
                if array is None:
                    dtype = 'float'
                elif hasattr(array,'_dtype'):
                    dtype = array._dtype
                elif hasattr(array,'dtype'):
                    dtype = fixType(array.dtype)
                else:
                    a = array[0]
                    while hasattr(a,'__len__'):
                        a = a[0]
                    if hasattr(array,'_dtype'):
                        dtype = a._dtype
                    elif hasattr(array,'dtype'):
                        dtype = fixType(a.dtype)
                    else:
                        dtype = fixType(type(a))
            self._dtype = dtype
            if type(array) is Epetra.MultiVector:
                if map is None:
                    self._mV = array
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
                    vLength = array.GlobalLength()
                    self._vLength = vLength
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
                    oldMap = array.Map()
                    DistToPers = Epetra.Import(map,oldMap)
                    mv = Epetra.MultiVector(map,array.NumVectors())
                    mv.Import(array, DistToPers, Epetra.Insert)
                    self._mV = mv
                    return
            if array is None:
                isTril = True
                if shape is None:
                    if map is None:
                        raise TypeError("_TrilinosArray.__init__() needs a map, shape, or array")
                    else:
                        shape = (map.NumGlobalElements(),)
                if vLength is None:
                    vLength = shape[-1]
                if map is None:
                    comm = Epetra.PyComm()
                    map = Epetra.Map(vLength,0,comm)
                else:
                    comm = map.Comm()
                    vLength = map.NumGlobalElements()
                shape = shape[:-1]+(vLength,)
                size = 1
                for i in shape[:-1]:
                    size *= i
                mv = Epetra.MultiVector(map,size)
            else:
                isTril = isinstance(array,_TrilinosArray)
                if not isTril:
                    narray = numpy.array(array,dtype='object')
                    t1 = narray.take([0])[0]
                    isTril = isinstance(t1,_TrilinosArray)
                    if isTril:
                        nshape = narray.shape
                        tshape = t1.shape
                        shape = nshape+tshape
                        tsize = 1
                        for k in tshape[:-1]:
                            tsize*=k
                        t1 = t1.multiVector
                        comm = t1.Comm()
                        vLength = t1.GlobalLength()
                        if map is None:
                            map = Epetra.Map(vLength,0,comm)
                        mv = Epetra.MultiVector(map,narray.size*tsize)
                        f = narray.flat
                        curN = 0
                        for i in f:
                            for k in i.multiVector:
                                mv[curN] = k
                                curN+=1
                else:
                    if shape is None or shape == array.shape:
                        shape = array.shape
                        array = array.multiVector
                        vLength = array.GlobalLength()
                        comm = array.Comm()
                        map = array.Map()
                        mv = Epetra.MultiVector(Epetra.Copy,array)
                    else:
                        array = array.multiVector
                        PersonalMap = Epetra.Map(-1, range(0, array.GlobalLength()), 0, v.Comm())
                        DistToPers = Epetra.Import(PersonalMap, array.Map())
                        PersonalV = Epetra.Vector(PersonalMap)
                        PersonalV.Import(array, DistToPers, Epetra.Insert)
                        vLength = shape[-1]
                        PersonalV = PersonalV.reshape(-1,vLength)
                        comm = array.comm
                        if map is None:
                            map = Epetra.Map(vLength,0,comm)
                        else:
                            vLength = map.NumGlobalElements()
                            shape = shape[:-1]+(vLength,)
                        if type(shape) != tuple:
                            shape = (shape,)
                        PersonalV = PersonalV[...,map.MyGlobalElements()]
                        mv = Epetra.MultiVector(map,narray)
            if not isTril:
                del t1
                if shape is None:
                    if array is not None:
                        shape = narray.shape
                    else:
                        shape = (map.NumGlobalElements(),)
                if vLength is None:
                    vLength = shape[-1]
                if map is None:
                    comm = Epetra.PyComm()
                    nvLength = vLength
                    map = Epetra.Map(vLength,0,comm)
                else:
                    comm = map.Comm()
                    nvLength = map.NumGlobalElements()
                    shape = shape[:-1]+(vLength,)
                try:
                    narray = narray.reshape(-1,vLength).astype('float')
                    narray = narray[...,map.MyGlobalElements()]
                except ValueError:
                    narray = narray.reshape(-1,narray.shape[-1]).astype('float')
                mv = Epetra.MultiVector(map,narray)
            self._shape = shape
            size = 1
            strides = ()
            for s in shape:
                strides = (size,)+strides
                size *= s
            self._strides = strides
            self._size = size
            self._vLength = vLength
            self._mV = mv
            self._dims = len(shape)
            inds = numpy.arange(size/shape[-1])
            if self._dims != 1:
                self._indices = inds.reshape(shape[:-1])
            else:
                self._indices = inds
            self._reprMV = numpy.arange(vLength)
    
    def take(self,ids,axis=None,mode=None):
        if debug:
            print "IN take\n"
        currMap = self.map
        comm = self.comm
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
                    resMap = ids.map
                    resVec = Epetra.MultiVector(resMap,1)
                    throughMap = Epetra.Map(-1,ids.array.astype(int),0,comm) #MUST be changed on addition of intvectors
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
                axis += dims
            if axis == dims-1:
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
        if debug:
            print "IN fillWith\n"
        self._mV.PutScalar(i)
    
    def copy(self):
        if debug:
            print "IN copy\n"
        newMV = Epetra.MultiVector(self.map,self.array.copy())
        newTril = _TrilinosArray(init=False)
        newTril._mV = newMV
        newTril._shape = self._shape
        newTril._vLength = self._vLength
        newTril._strides = self._strides
        newTril._size = self._size
        newTril._dims = self._dims
        newTril._indices = self._indices
        newTril._reprMV = self._reprMV
        newTril._dtype = self._dtype
        newTril.takes = self.takes
        return newTril
    
    def setMap(self,m):
        if debug:
            print "IN setMap\n"
        newMV = Epetra.MultiVector(m,self._indices.size)
        imp = Epetra.Import(m,self.map)
        newMV.Import(self._mV,imp,Epetra.Insert)
        self._vLength = m.NumGlobalElements()
        self._reprMV = numpy.arange(self._vLength)
        self._shape = self._shape[:-1]+(self._vLength,)
        self._mV = newMV
    
    def reshape(self,shape,*args,**copy):
        if debug:
            print "IN reshape\n"
        if copy == {}:
            copy['copy'] = True
        elif not copy.has_key('copy'):
            raise TypeError, "Invalid keyword arguments.  The only valid keyword is 'copy'."
        shape = (args and (shape,)+args) or shape
        if copy['copy']:
            newTril = _TrilinosArray(init=False)
            newTril._mV = self._mV
            newTril._shape = shape
            newTril._vLength = self._vLength
            size = 1
            strides = ()
            for s in shape:
                strides = (size,)+strides
                size *= s
            newTril._strides = strides
            newTril._size = self._size
            newTril._dims = len(shape)
            newTril._indices = self._indices.reshape(shape[:-1])
            newTril._reprMV = self._reprMV
            newTril._dtype = self._dtype
            return newTril
        else:
            self._shape = shape
            self._indices = self._indices.reshape(shape[:-1])
            self._dims = len(shape)
    
    def astype(self,dtype):
        if debug:
            print "IN astype\n"
        res = self.copy()
        dtype = fixType(dtype)
        res._dtype=dtype
        return res
    
    def __getslice__(self,i,j):
        if debug:
            print "IN __getslice__\n"
        return self.__getitem__((slice(i,j,None),))
    
    def __getitem__(self,y):
        if debug:
            print "IN __getitem__"
            print "self = "+repr(self)+"y = "+str(y)+"\n"
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
        m = self.map
        myInds = m.MyGlobalElements()
        comm = self.comm
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
                    return res.astype(self._dtype)
                elif type(forMV) is _TrilinosArray:
                    print repr(forMV)
                    els = mv[0,forMV.array]
                    numEls = els.size
                    totEls = comm.SumAll(numEls)
                    m = Epetra.Map(-1,numEls,0,comm)
                    newMV = Epetra.MultiVector(m,els)
                    return _TrilinosArray(newMV,shape=totEls)
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
                        return res.astype(self._dtype)
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
                if forMV == Ellipsis or (type(forMV) == slice and forMV.stop >= vLength and forMV.start <= 0):
                    inds = indices[y]
                    newMV = mv[inds.reshape(-1)]
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
        if debug:
            print "IN __setslice__\n"
        self.__setitem__((slice(i,j,None),),a)
    
    def __setitem__(self,y,a):
        if debug:
            print "IN __setitem__\n"
        t = False
        if type(a) == _TrilinosArray:
            a = a.__array__()
            t = True
        if type(y) == _TrilinosArray:
            y = y.__array__()
            t = True
        if t:
            self._mV[y] = a
        if type(y) is not tuple:
            y = (y,)
        print repr(y),repr(a)
        m = self.map
        comm = self.comm
        dims = self._dims
        indices = self._indices
        test = self._reprMV
        myInds = self.map.MyGlobalElements()
        mv = self._mV
        if len(y) == dims or y.__contains__(Ellipsis):
            indices = indices[y[:-1]]
            test = test[y[-1]]
            if type(test) is not numpy.ndarray:
                test = numpy.array([test])
            test = [m.LID(i) for i in test if i in myInds]
        else:
            indices = indices[y]
        if dims is 1:
            if hasattr(a,'shape') and len(a.shape) is 2:
                 a = a[0]
            mv[0,test] = a
        else:
            mv[indices][:,test] = a
            
    def __mul__(self,other):
        if debug:
            print "IN __mul__\n"
        if type(other) == _TrilinosArray:
            dtype = other._dtype
        elif type(other) == numpy.ndarray:
            dtype = fixType(other.dtype)
        else:
            dtype = fixType(type(other))
        res = self.astype(retType(self._dtype,dtype))
        if hasattr(other,'_mV'):
            res._mV[:] *= other._mV
            return res
        res._mV[:] *= other
        return res
    
    def __add__(self,other):
        if debug:
            print "IN __add__\n"
        if type(other) == _TrilinosArray:
            dtype = other._dtype
        elif type(other) == numpy.ndarray:
            dtype = fixType(other.dtype)
        else:
            dtype = fixType(type(other))
        res = self.astype(retType(self._dtype,dtype))
        if hasattr(other,'_mV'):
            res._mV[:] += other._mV
            return res
        res._mV[:] += other
        return res

    def __div__(self,other):
        if debug:
            print "IN __div__\n"
        if type(other) == _TrilinosArray:
            dtype = other._dtype
        elif type(other) == numpy.ndarray:
            dtype = fixType(other.dtype)
        else:
            dtype = fixType(type(other))
        res = self.astype(retType(self._dtype,dtype))
        if hasattr(other,'_mV'):
            res._mV[:] /= other._mV
            return res
        res._mV[:] /= other
        return res

    def __sub__(self,other):
        if debug:
            print "IN __sub__\n"
        if type(other) == _TrilinosArray:
            dtype = other._dtype
        elif type(other) == numpy.ndarray:
            dtype = fixType(other.dtype)
        else:
            dtype = fixType(type(other))
        res = self.astype(retType(self._dtype,dtype))
        if hasattr(other,'_mV'):
            res._mV[:] -= other._mV
            return res
        res._mV[:] -= other
        return res

    def __rmul__(self,other):
        if debug:
            print "IN __rmul__\n"
        if type(other) == _TrilinosArray:
            dtype = other._dtype
        elif type(other) == numpy.ndarray:
            dtype = fixType(other.dtype)
        else:
            dtype = fixType(type(other))
        res = self.astype(retType(self._dtype,dtype))
        if hasattr(other,'_mV'):
            res._mV[:] *= other._mV
            return res
        res._mV[:] *= other
        return res
    
    def __radd__(self,other):
        if debug:
            print "IN __radd__\n"
        if type(other) == _TrilinosArray:
            dtype = other._dtype
        elif type(other) == numpy.ndarray:
            dtype = fixType(other.dtype)
        else:
            dtype = fixType(type(other))
        res = self.astype(retType(self._dtype,dtype))
        if hasattr(other,'_mV'):
            res._mV[:] += other._mV
            return res
        res._mV[:] += other
        res._dtype = retType(self._dtype,dtype)
        return res
    
    def __rdiv__(self,other):
        if debug:
            print "IN __rdib__\n"
        if hasattr(other,'_mV'):
            res = other.astype(retType(self._dtype,other._dtype))
            res._mV[:] /= self._mV
            return res
        res = self.__inv__()
        res._mV[:] *= other
        if type(other) == numpy.ndarray:
            dtype = other.dtype
        else:
            dtype = type(other)
        dtype = fixType(dtype)
        res._dtype = retType(self._dtype,dtype)
        return res

    def __rsub__(self,other):
        if debug:
            print "IN __rsub__\n"
        if hasattr(other,'_mV'):
            res = other.astype(retType(self._dtype,other._dtype))
            res._mV[:] -= self._mV
            return res
        res = self.__neg__()
        res._mV[:] += other
        if type(other) == numpy.ndarray:
            dtype = other.dtype
        else:
            dtype = type(other)
        dtype = fixType(dtype)
        res._dtype = retType(self._dtype,dtype)
        return res

    def __imul__(self,other):
        if debug:
            print "IN __imul__\n"
        if hasattr(other,'_mV'):
            self._mV[:] *= other._mV
            self._dtype = retType(self._dtype,other._dtype)
            return self
        self._mV[:] *= other
        if type(other) == numpy.ndarray:
            dtype = other.dtype
        else:
            dtype = type(other)
        dtype = fixType(dtype)
        self._dtype = retType(self._dtype,dtype)
        return self
    
    def __iadd__(self,other):
        if debug:
            print "IN __iadd__\n"
        if hasattr(other,'_mV'):
            self._mV[:] += other._mV
            self._dtype = retType(self._dtype,other._dtype)
            return self
        self._mV[:] += other
        if type(other) == numpy.ndarray:
            dtype = other.dtype
        else:
            dtype = type(other)
        dtype = fixType(dtype)
        self._dtype = retType(self._dtype,dtype)
        return self
    
    def __idiv__(self,other):
        if debug:
            print "IN __idib__\n"
        if hasattr(other,'_mV'):
            self._mV[:] /= other._mV
            self._dtype = retType(self._dtype,other._dtype)
            return self
        self._mV[:] /= other
        if type(other) == numpy.ndarray:
            dtype = other.dtype
        else:
            dtype = type(other)
        dtype = fixType(dtype)
        self._dtype = retType(self._dtype,dtype)
        return self

    def __isub__(self,other):
        if debug:
            print "IN __isub__\n"
        if hasattr(other,'_mV'):
            self._mV[:] -= other._mV
            self._dtype = retType(self._dtype,other._dtype)
            return self
        self._mV[:] -= other
        if type(other) == numpy.ndarray:
            dtype = other.dtype
        else:
            dtype = type(other)
        dtype = fixType(dtype)
        self._dtype = retType(self._dtype,dtype)
        return self
    
    def __inv__(self):
        if debug:
            print "IN __inv__\n"
        res = self.astype('float')
        res._mV[:] = 1/res._mV
        return res

    def __neg__(self):
        if debug:
            print "IN __neg__\n"
        res = self.astype(retType(self._dtype,'int'))
        res._mV[:] = -res._mV
        return res
    
    def __or__(self,other):
        if debug:
            print "IN __or__\n"
        a = self.__array__()
        if hasattr(other,'_mV'):
            other = other.__array__()
        a |= other
        mv = Epetra.MultiVector(self.map,a)
        return _TrilinosArray(mv,shape=self._shape,dtype = fixType(a.dtype))
    
    def __eq__(self,other):
        if debug:
            print "IN __eq__\n"
        a = self.__array__()
        if hasattr(other,'_mV'):
            other = other.__array__()
        a = a == other
        mv = Epetra.MultiVector(self.map,a)
        return _TrilinosArray(mv,shape=self._shape,dtype=fixType(a.dtype))
    
    def __iter__(self):
        if debug:
            print "IN __iter__\n"
        if self._dims is 1:
            return self._mV.array[0].__iter__()
        return self._mV.__iter__()
    
    def __len__(self):
        if debug:
            print "IN __len__\n"
        if self._dims is 1:
            return self.vLength
        return self._mV.__len__()

    def __array__(self,dtype=None):
        if dtype is None:
            dtype = self._dtype
        else:
            dtype = fixType(dtype)
        if debug:
            print "IN __array__\n"
        if self._dims is 1:
            return self._mV.array[0].astype(dtype)
        return self._mV.array.astype(dtype)
    
    def __array_wrap__(self,obj):
        if debug:
            print "IN __array_wrap__\n"
        return _TrilinosArray(obj,shape = self._shape,map = self.map,dtype = self._dtype)
    
    def __str__(self):
        if debug:
            print "IN __str__\n"
        return self.array.reshape(self.shape[:-1]+(-1,)).__str__()

    def __repr__(self):
        if debug:
            print "IN __repr__\n"
        return "_TrilinosArray("+self.multiVector.__str__()+" dtype = "+self._dtype+" shape = "+self.shape.__str__()+")"

    if debug:
        def __getattribute__(self,attr):
            print attr
            return object.__getattribute__(self,attr)

    shape = property(lambda self:self._shape[:-1]+(self._vLength,),reshape,doc="The global shape of the array")
    vLength = property(lambda self:self._vLength,doc="The global length of each vector in the MultiVector")
    map = property(lambda self:self._mV.Map(),setMap,doc="The Epetra map used to distribute the MultiVector")
    comm = property(lambda self:self._mV.Comm(),doc="The Epetra communicator used by this _TrilinosArray")
    multiVector = property(lambda self:self._mV,doc="The multiVector wrapped by this _TrilinosArray")
    array = property(__array__,doc="The array representing the local _TrilinosArray")
    dtype = property(lambda self:numpy.dtype(self._dtype))

    #def __getattr__(self,attr):
    #    print attr
    #    mv = object.__getattribute__(self,'_mV')
    #    attribute = getattr(mv,attr)
    #    if callable(attribute):
    #        shape = object.__getattribute__(self,'_shape')
    #        map = mv.Map()
    #        return lambda *args,**kwargs: _wrap(attribute,shape,map,*args,**kwargs)
    #    else:
    #        return attribute

def isTrilArray(array):
    return type(array)==_TrilinosArray

def _wrap(fn,retShape,retmap,*args,**kwargs):
    obj = fn(*args,**kwargs)
    if type(obj) == Epetra.MultiVector:
        return _TrilinosArray(obj,shape=retShape,map=retmap)
    return obj

def arange(stop,start = None,step = 1,shape=None,map = None,dtype='int'):
    if map is None:
        if shape is None:
            if step is None:
                if start is None:
                    return _TrilinosArray(numpy.arange(stop))
                return _TrilinosArray(numpy.arange(stop,start))
            return _TrilinosArray(numpy.arange(stop,start,step))
        if step is None:
            if start is None:
                return _TrilinosArray(numpy.arange(stop).reshape(shape[:-1]+(-1,)))
            return _TrilinosArray(numpy.arange(stop,start).reshape(shape[:-1]+(-1,)))
        return _TrilinosArray(numpy.arange(stop,start,step).reshape(shape[:-1]+(-1,)))
    if start is None:
        stop = start
        start = 0
    comm = map.Comm()
    nme = map.NumMyElements()
    myStop = comm.ScanSum(nme)
    myStart = myStop-nme
    arr = numpy.arange(start+myStart*step,start+myStop*step,step,dtype='d')
    mv = Epetra.MultiVector(map,arr)
    if shape is not None:
        return _TrilinosArray(mv,shape=shape,dtype=dtype)
    else:
        return _TrilinosArray(mv,shape=map.NumGlobalElements(),dtype=dtype)

def fixType(t):
    if str(t).count('int')>0 or t is int \
           or t in ['l','i']:
        return 'int'
    if str(t).count('float')>0 or str(t).count('double')>0 \
           or t is float or t in ['f','d']:
        return 'float'
    if str(t).count('bool')>0 or t is bool \
           or t in ['b1']:
        return 'bool'
    return t

def retType(type1,type2):
    if type1 == 'float' or type2 == 'float':
        return 'float'
    if type1 == 'int' or type2 == 'int':
        return 'int'
    return 'bool'

