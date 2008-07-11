#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "trilinosArray.py"
 #                                    created: 7/3/08 {10:23:17 AM} 
 #                                last update: 7/10/08 {11:45:26 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Olivia Buzek   <olivia.buzek@nist.gov>
 #  Author: Daniel Stiles  <dastiles@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##
 
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
    def __init__(self, array=None, shape=None, map=None, dType='l', \
                 parallel=True):
        """
        Creates a trilArr

        :Parameters:
          - `shape`:the shape of the array.  If passed with an array, it overrides the shape of the array.
          - `map`:an Epetra.Map or Epetra.BlockMap describing how to split the job between processors
          - `dType`:the type of data in the array.  However, due to Trilinos limitations, double will be converted to float, and anything besides float or double will become long
          - `parallel`:whether or not this array should be parallelized
          - `array`:the array to be used.  Any iterable type is accepted here (numpy.array, list, tuple).  Default is all zeros.  The data type of the array overrides dType if it is passed in.
        """
        import operator
        if shape is None and array is None:
            print "FAIL: Must specify either shape or vector."

        if array is None or str(type(array)).count("Epetra") == 0:
            
            if array is not None:
                self.shp  = trilShape(numpy.array(array).shape)
                if str(numpy.array(array).dtype).count("float") != 0:
                    dType = 'f'
            if shape is not None:
                self.shp = trilShape(shape)

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
                        tmpArray = numpy.array(array).reshape(-1)
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
                    tmpArray = numpy.array(array).reshape(-1)
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
            if str(numpy.array(array).dtype).count("float") != 0:
                dType = 'f'
            self.vector = array
            self.comm = array.Comm()
            self.eMap = array.Map()
            self.shp = trilShape(self.eMap.NumGlobalElements())
            self.shp.setMap(self.eMap)
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
          - `ids`: Where to place the values
          - `values`: The values to insert.  If there are less than there are ids, loops through the list multiple times

            >>> t = trilArr(shape=(4,))
            >>> t.insertValues([0],[5])
            >>> t.allElems()
            trilArr([5, 0, 0, 0])
            >>> t.insertValues([1,2,3],[7,8])
            >>> t.allElems()
            trilArr([5, 7, 8, 7])
        """

        if self.eMap is not None:
            elms = list(self.eMap.MyGlobalElements())
            if type(values) != int:
                values = [v for (i,v) in zip(ids,values*((len(ids)+1)/len(values))) if elms.count(i)>0]
            ids = [self.eMap.LID(i) for i in ids if list(elms).count(i)>0]
        numpy.put(self.array, ids, values)
    
    def take(self,ids):
        """
        Takes values out of the array
        
        :Parameters:
          - `ids`: What values to take
        
            >>> t = trilArr(array=[1,2,3,4])
            >>> t.take([1,2])
            array([2, 3])
         """
        return self.globalTake(ids)

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
        allEls = [l for (el,proc) in zip(allEls,range(procs)) \
                  for (l,pos) in zip(el,range(sizes[proc]))]
        allEls = numpy.array(allEls).reshape(shape)
        return allEls

    def localTake(self, ids):
        indices = numpy.array(ids)
        indices = indices.reshape(-1)
        glob = self.eMap.MyGlobalElements()
        num = self.eMap.NumGlobalElements()
        for (ind,i) in zip(indices,range(len(indices))):
            if ind<0:
                indices[i] = ind+num
        myIDs = [self.eMap.LID(el) for el in indices \
                 if list(glob).count(el)>=1]
        if myIDs == []: return []
        return self.vector[myIDs]

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
            >>> t.arccos().allElems()
            trilArr([ 0.,  0.,  0.,  0.])
        """
        return self._applyFloatFunction(numpy.arccos)

    def arccosh(self):
        """
        arccosh of this array

            >>> t = trilArr(array=[1, 1, 1, 1])
            >>> t.arccosh().allElems()
            trilArr([ 0.,  0.,  0.,  0.])
        """
        return self._applyFloatFunction(numpy.arccosh)

    def arcsin(self):
        """
        arccos of this array

            >>> t = trilArr(array=[1, 1, 1, 1])
            >>> t.arcsin().allElems()
            trilArr([ 1.57079633,  1.57079633,  1.57079633,  1.57079633])
        """
        return self._applyFloatFunction(numpy.arcsin)

    def arcsinh(self):
        """
        arcsinh of this array

            >>> t = trilArr(array=[1, 1, 1, 1])
            >>> t.arcsinh().allElems()
            trilArr([ 0.88137359,  0.88137359,  0.88137359,  0.88137359])
        """
        return self._applyFloatFunction(numpy.arcsinh)

    def arctan(self):
        """
        arctan of this array

            >>> t = trilArr(array=[1, 1, 1, 1])
            >>> t.arctan().allElems()
            trilArr([ 0.78539816,  0.78539816,  0.78539816,  0.78539816])
        """
        return self._applyFloatFunction(numpy.arctan)

    def arctanh(self):
        """
        arctanh of this array

            >>> t = trilArr(array=[.5, .5, .5, .5])
            >>> t.arctanh().allElems()
            trilArr([ 0.54930614,  0.54930614,  0.54930614,  0.54930614])
        """
        return self._applyFloatFunction(numpy.arctanh)

    def arctan2(self, other):
        """
        arctan of this array/other

        :Parameters:
          - `other`: The array in the denominator

            >>> n = trilArr(array=[0, 0, 0, 0])
            >>> d = trilArr(array=[1, 1, 1, 1])
            >>> n.arctan2(d).allElems()
            trilArr([ 0.,  0.,  0.,  0.])
            >>> d.arctan2(n).allElems()
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

    def reshape(self, *args, **kwargs):
        ## reshape checks need to be done
        ## before a copy is made
        if kwargs.has_key("copy"):
            copy = kwars["copy"]
        else:
            if len(args) > 0 and type(args[-1]) == bool:
                copy = args[-1]
                args = args[:-1]
            else:
                copy = False
        if kwargs.has_key("shape"):
            copy = kwars["shape"]
        else:
            if len(args) > 0:
                if len(args) == 1 and type(args[0]) == tuple:
                    shape = args[0]
                else:
                    shape = args
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
        return self.__getitem__(slice(i,j,None))
    
    def __setitem__(self, i, y):
        i = self.shp.getLocalIndex(i)
        self.vector.__setitem__(i[0], y)

    def __getitem__(self, y):
        y = self.shp.getLocalIndex(y)
        a = self.vector.__getitem__(y[0])
        s = y[1]
        if s == (): return a[0]
        return trilArr(array = a,shape = s)

    def __copy__(self):
        if self.eMap.NumMyElements() == self.eMap.NumGlobalElements():
            plell = false
        else:
            plell = true
        return trilArr(array = self.vector.copy(), shape = self.shp.copy(), \
                       dType = self.dType, parallel = plell)

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

    def __len__(self):
        return self.shp.globalShape[0]

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
        tmp.reverse()
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
        return (self._globalToLocal(ind[0]),ind[1])

    def _globalToLocal(self, i):
        if self.map is None:
            return -1
        if type(i)==int:
            return self.map.LID(i)
        else:
            return [self.map.LID(j) for j in i]

    def _intToSlice(self, i):
        if type(i)==slice:
            return i
        elif i == -1:
            return slice(i,None,None)
        else:
            return slice(i,i+1,None)

    def _fillToDim(self, i):
        i = list(i)
        while len(i)<self.dimensions:
            i.append(slice(None,None,None))
        return tuple(i)
    
    def _globalTranslateIndices(self, index):
        if type(index)==int or type(index)==list:
            index=[(index,)]
        elif type(index)==slice:
            index=[(index,)]
        elif type(index)==tuple:
            if type(index[0])!=int and type(index[0])!=slice and type(index[0])!=list:
                while type(index)!=int and len(index)==1:
                    index=index[0]
                if type(index)==int:
                    index=[(index,)]
                elif len(index)<=self.dimensions:
                    index = [[i[el] for i in index] for el in range(len(i))]
            else:
                index = [tuple(index)]
        index = [self._fillToDim(i) for i in index]
        index = [el for i in index for el in self._globalTranslateSlices(i)]
        s = index[1]
        index = index[0]
        indices = []

        for ind in index:
            if self._dimensions(ind)>self.dimensions:
                return -1
            if not sum([i<j for (i,j) in zip(ind,self.globalShape)]):
                return -2

            lineIndex = 0
            for (mult,i) in zip(self.steps,range(len(ind))):
                lineIndex += mult*ind[i]

            indices.append(lineIndex)

        return (indices,s)

    def _globalTranslateSlices(self, sls):
        sls = list(sls)
        o = list(sls)
        for (el,i) in zip(sls,range(len(sls))):
            if type(el)==int:
                sls[i]=self._intToSlice(el)
        dims = [numpy.arange(i) for i in self.globalShape]
        res = [tuple(list(dim[sl])) for (sl,dim) in zip(sls,dims)]
        s = tuple([len(d) for (d,n) in zip(res,o) if str(type(n)).count("int")==0])
        k = [len(i) for i in res]
        k2 = [len(i) for i in res]
        for i in range(len(k2))[1:]:
            k2[i]*=k2[i-1]
            if k2[i] == 0:
                return [()]

        m = k2[-1]

        ans = [tuple([p for el in \
                      [(tup[i],)*(m/l) for i in range(j)]\
                      for p in el]) \
               for (tup,l,j) in zip(res,k2,k)]

        k2 = [1] + k2

        fin = [tuple([p for el in (tup,)*z \
                      for p in el]) \
               for (tup,z) in zip(ans,k2[:-1])]

        inds = [tuple([i[j] for i in fin]) for j in range(m)]
        return (inds,s)

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
        return shape

    def reshape(self, shape):
        shape = self._shapeCheck(shape)
        if self._size(shape) < 0:
            un = -1
            tot = 1
            for i in shape:
                if i == -1:
                    if un < 0:
                        un = i
                    else:
                        print "ERROR: Only one unspecified dimension is allowed."
                        return -1
                elif i > 0:
                    tot *= i
                else:
                    print "ERROR: Negative sizes are not allowed."
                    return -1
            p = self.actualShape*1./tot
            if numpy.ceil(p) != numpy.floor(p):
                print "ERROR: New shape doesn't fit"
                return -1
            shape = list(shape)
            shape[un] = int(p)
            shape = tuple(shape)

        if self.actualShape != self._size(shape):
            print "ERROR: New shape is differently sized from old shape."
            return -1
        self.globalShape = shape
        self.actualShape = self._size(shape)
        self.dimensions = self._dimensions(shape)

        mult = 1
        tmp = []
        for i in range(len(self.globalShape)+1)[1:]:
            tmp.append(mult)
            mult *= self.globalShape[-i]
        tmp.reverse()
        self.steps = tuple(tmp)

        return 1

    def __str__(self):
        return self.globalShape.__str__()

    def __repr__(self):
        return "trilShape("+self.globalShape.__repr__()+")"

    def __copy__(self):
        return trilShape(self.globalShape, self.eMap)
    

def isTrilArray(obj):
    return isinstance(obj, trilArr)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
