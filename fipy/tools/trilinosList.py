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
                epetraMap = Epetra.Map(vLength,0,comm)
                comm = Epetra.PyComm()
            else:
                comm = epetraMap.Comm()
                vLength = epetraMap.NumGlobalElements()
                shape = shape[:-1]+(vLength,)
            self._shape = shape
            self._vLength = vLength
            self._comm = comm
            self._map = epetraMap
            self._mV = Epetra.MultiVector(epetraMap)
        else:
            flag = isinstance(array,_TrilinosArray) or type(array) == Epetra.MultiVector
            if not flag:
                narray = numpy.array(array,dtype='object')
                t1 = narray.take([0])[0]
                flag = flag or isinstance(t1,_TrilinosArray) or type(t1) == Epetra.MultiVector
                if flag:
                    depth = len(narray.shape)
                    shape = narray.shape + t1.shape
                    
            else:
                if shape is None or shape == array.shape or shape == (array.shape,):
                    self._shape = array.shape
                    if isinstance(array,_TrilinosArray):
                        array = array.multiVector
                    self._vLength = array.GlobalLength()
                    comm = array.Comm()
                    self._comm = comm
                    oldMap = array.Map()
                    if epetraMap is None:
                        epetraMap = oldMap
                    self._map = epetraMap
                    DistToPers = Epetra.Import(epetraMap,oldMap)
                    PersonalV = Epetra.MultiVector(epetraMap,array.NumVectors())
                    PersonalV.Import(array, DistToPers, Epetra.Insert)
                    self._mV = PersonalV
                else:
                    if isinstance(array,_TrilinosArray):
                        array = array.multiVector
                    PersonalMap = Epetra.Map(-1, range(0, array.GlobalLength()), 0, v.Comm())
                    DistToPers = Epetra.Import(PersonalMap, array.Map())
                    PersonalV = Epetra.Vector(PersonalMap)
                    PersonalV.Import(array, DistToPers, Epetra.Insert)
                    vLength = shape[-1]
                    PersonalV = PersonalV.reshape(-1,vLength)
                    self._comm = array.comm
                    if epetraMap is None:
                        epetraMap = Epetra.Map(vLength,0,comm)
                    else:
                        vLength = epetraMap.NumGlobalElements()
                        shape = shape[:-1]+(vLength,)
                    if type(shape) != tuple:
                        shape = (shape,)
                    self._shape = shape
                    self._vLength = vLength
                    self._map = epetraMap
                    PersonalV = PersonalV[...,epetraMap.MyGlobalElements()]
                    self._mV = Epetra.MultiVector(epetraMap,narray)
                    
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
                epetraMap = Epetra.Map(vLength,0,comm)
                comm = Epetra.PyComm()
            else:
                comm = epetraMap.Comm()
                vLength = epetraMap.NumGlobalElements()
                shape = shape[:-1]+(vLength,)
            self._shape = shape
            self._vLength = vLength
            self._comm = comm
            self._map = epetraMap
            narray = narray[...,epetraMap.MyGlobalElements()]
            self._mV = Epetra.MultiVector(epetraMap,narray)
        
    def __str__(self):
        return self.multiVector.__str__()

    def __repr__(self):
        return "TrilinosArray("+self.multiVector.__str__()+" shape = "+self.shape.__str__()+")"

class trilIntArr:
    def __init__(self,array=None,shape=None,map=None,dtype=None):
        pass


def __call__(*args,**kwargs):
    return _TrilinosArray(args,kwargs)
