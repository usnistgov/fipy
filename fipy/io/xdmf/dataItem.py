#!/usr/bin/env python
## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "dataItem.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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

import os
from xml.dom import minidom

from fipy.tools import numerix

from fipy.io.xdmf import _xml_to_array
from fipy.io.xdmf.node import _Node

class DataItem(_Node):
    @classmethod
    def from_xml(cls, document, node):
        format = node.getAttribute("Format") or "XML"
        if format == "XML":
            data = XMLDataItem(document=document, node=node)
        elif format == "HDF":
            data = HDFDataItem(document=document, node=node)
        else:
            raise Exception("Unknown DataItem Format: ", format)
            
        return data
        
    @staticmethod
    def _node(document, arr, name=None):
        arr = numerix.asarray(arr)
        data = document.node.createElement("DataItem")
        if name is not None:
            data.setAttribute("Name", name)
        if arr.shape == ():
            data.setAttribute("Dimensions", "1")
        else:
            data.setAttribute("Dimensions", " ".join(str(v) for v in arr.shape))
            
        return data

class XMLDataItem(DataItem):
    @classmethod
    def from_array(cls, document, arr, name=None):
        data = cls._node(document=document, arr=arr, name=name)
        data.setAttribute("Format", "XML")
                                   
        arr = numerix.asarray(arr)
        if arr.shape == ():
            text = document.node.createTextNode(str(arr))
        else:
            text = document.node.createTextNode(" ".join(str(v) for v in arr.flat))
        data.appendChild(text)
        
        return cls(document=document, node=data)
        
    @property
    def array(self):
        self.node.normalize()
        
        shape = _xml_to_array(self.node.getAttribute("Dimensions"))
        return _xml_to_array(self.node.firstChild.data).reshape(shape)
        
class HDFDataItem(DataItem):
    @classmethod
    def from_array(cls, document, arr, name, h5filename):
        data = cls._node(document=document, arr=arr, name=name)
        data.setAttribute("Format", "HDF")
                                   
        arr = numerix.asarray(arr)
        
        import h5py
        f = h5py.File(h5filename, mode='a')
        if name in f:
            del f[name]
        f.create_dataset(name, data=arr)
        f.close()
        
        text = document.node.createTextNode(os.path.split(h5filename)[1] + ":/" + name)
        data.appendChild(text)
        
        return cls(document=document, node=data)

    @property
    def array(self):
        self.node.normalize()
        h5filename, name = self.node.firstChild.data.split(':')
        
        import h5py
        f = h5py.File(h5filename.strip(), mode='r')
        arr = numerix.asarray(f[name.strip()])
        f.close()
        
        return arr
