#!/usr/bin/env python
## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "document.py"
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

from fipy.io.xdmf.dataItem import DataItem, HDFDataItem, XMLDataItem
from fipy.io.xdmf.node import _Node
from fipy.variables.variable import Variable
from fipy.variables.cellVariable import CellVariable

class Attribute(_Node):
    @classmethod
    def from_Variable(cls, document, var):
        attribute = document.node.createElement("Attribute")
        attribute.setAttribute("Name", var.name)
        attribute.setAttribute("AttributeType", "Scalar")
        attribute.setAttribute("Center", "Grid")

        data = XMLDataItem.from_array(document=document, arr=var.value)
        attribute.appendChild(data.node)

        return cls(document=document, node=attribute)
        
    @property
    def data(self):
        data = self.node.getElementsByTagName("DataItem")[0]
        return DataItem.from_xml(document=self.document, node=data)
    
    @property
    def variable(self):
        name = self.node.getAttribute("Name")
        return Variable(name=name, value=self.data.array)

class MeshAttribute(Attribute):
    def __init__(self, document, node, grid=None):
        Attribute.__init__(self, document=document, node=node)
        self.grid = grid
        
    @classmethod
    def from_MeshVariable(cls, document, var, grid, h5filename):
        from fipy.variables import CellVariable, FaceVariable
        if isinstance(var, CellVariable):
            attr = CellAttribute.from_CellVariable(document=document, var=var, grid=grid, h5filename=h5filename)
        elif isinstance(var, FaceVariable):
            attr = FaceAttribute.from_FaceVariable(document=document, var=var, grid=grid, h5filename=h5filename)
        else:
            raise Exception("Can't make Attribute from ", var.__class__.__name__)
            
        return attr
            
    @staticmethod
    def _node(document, var, data, h5filename):
        attribute = document.node.createElement("Attribute")
        attribute.setAttribute("Name", var.name)
        rank = ["Scalar", "Vector", "Tensor"]
        attribute.setAttribute("AttributeType", rank[var.rank])
        
        if data.size < 10:
            data = XMLDataItem.from_array(document=document, arr=data)
        else:
            data = HDFDataItem.from_array(document=document, arr=data, name=var.name, h5filename=h5filename)
        attribute.appendChild(data.node)
        
        return attribute

class CellAttribute(MeshAttribute):
    @classmethod
    def from_CellVariable(cls, document, var, grid, h5filename):
        data = grid.reshape_cells(var.value.copy())
        attribute = cls._node(document=document, var=var, data=data, h5filename=h5filename)
        attribute.setAttribute("Center", "Cell")

        return cls(document=document, node=attribute)
        
    @property
    def variable(self):
        name = self.node.getAttribute("Name")
        data = self.grid.unshape_cells(self.data.array)
        return CellVariable(mesh=self.grid.mesh, name=name, value=data)
        
class FaceAttribute(MeshAttribute):
    @classmethod
    def from_FaceVariable(cls, document, var, grid, h5filename):
        data = grid.reshape_faces(var.value.copy())
        attribute = cls._node(document=document, var=var, data=data, h5filename=h5filename)
        attribute.setAttribute("Center", "Face")

        return cls(document=document, node=attribute)
        
    @property
    def variable(self):
        name = self.node.getAttribute("Name")
        data = self.grid.unshape_faces(self.data.array)
        return FaceVariable(mesh=self.grid.mesh, name=name, value=data)
