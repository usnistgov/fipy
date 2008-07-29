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
        if array is None:
            flag = False
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
            flag = isinstance(array,_TrilinosArray) or type(array) == Epetra.MultiVector
            if not flag:
                narray = numpy.array(array,dtype='object')
                t1 = narray.take([0])[0]
                flag = isinstance(t1,_TrilinosArray) or type(t1) == Epetra.MultiVector
                if flag:
                    depth = len(narray.shape)
                    nshape = narray.shape
                    tshape = t1.shape
                    tsize = tshape[0]
                    for k in tshape[1:-1]:
                        tsize*=k
                    if isinstance(t1,_TrilinosArray):
                        t1 = t1.multiVector
                    comm = t1.Comm()
                    vLength = t1.GlobalLength()
                    if epetraMap is None:
                        epetraMap = Epetra.Map(vLength,0,comm)
                    mv = Epetra.MultiVector(epetraMap,narray.size*tsize)
                    f = narray.flat
                    curN = 0
                    for i in range(len(f)):
                        v = f[i]
                        if isinstance(v,_TrilinosArray):
                            v = v.multiVector
                        for k in range(len(v)):
                            mv[curN,:] = v[k,:]
                            curN+=1
            else:
                if shape is None or shape == array.shape or shape == (array.shape,):
                    shape = array.shape
                    if isinstance(array,_TrilinosArray):
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
                    if isinstance(array,_TrilinosArray):
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
                    
        if not flag:
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
    
    
    
    def __str__(self):
        return self.multiVector.__str__()

    def __repr__(self):
        return "_TrilinosArray("+self.multiVector.__str__()+" shape = "+self.shape.__str__()+")"

class trilIntArr:
    def __init__(self,array=None,shape=None,map=None,dtype=None):
        pass


def __call__(*args,**kwargs):
    return _TrilinosArray(args,kwargs)
