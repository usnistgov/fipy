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

from xml.dom import minidom

from fipy.tools import numerix

from fipy.io.xdmf.attribute import Attribute, CellAttribute
from fipy.io.xdmf.dataItem import DataItem, XMLDataItem, xml_to_array
from fipy.io.xdmf.node import _Node

class Grid(_Node):
    """Encapsulates a FiPy Mesh and the Variables declared on that Mesh
    """
    @classmethod
    def from_xml(cls, document, node):
        gridType = node.getAttribute("GridType")
        if gridType == "Collection":
            grid = CollectionGrid(document=document, node=node)
        else:
            topology = node.getElementsByTagName("Topology")[0]
            topologyType = topology.getAttribute("TopologyType") or topology.getAttribute("Type")
            if topologyType == "3DCoRectMesh":
                grid = ThreeDCoRectMeshGrid(document=document, node=node, mesh=None)
            else:
                raise Exception("Unknown TopologyType: " + topologyType)
            
        return grid
        
    def reshape_cells(self, data):
        return data
        
class MeshGrid(Grid):
    def __init__(self, document, node, mesh):
        _Node.__init__(self, document=document, node=node)
        self._mesh = mesh
        
    @property
    def variables(self):
        vars = []
        for node in self.node.childNodes:
            if node.nodeType == minidom.Node.ELEMENT_NODE:
                if node.tagName == "Attribute":
                    center = node.getAttribute("Center") or "Node"
                    if center == "Cell":
                        vars.append(CellAttribute(document=self.document, node=node, grid=self).variable)
                    elif center == "Face":
                        vars.append(FaceAttribute(document=self.document, node=node, grid=self).variable)
                    else:
                        raise Exception("Unknown Attribute center: ", center)

        return vars
   
class ThreeDCoRectMeshGrid(MeshGrid):
    def reshape_cells(self, data):
#         return data.reshape((-1,) + self.mesh.shape)
        data = numerix.rollaxis(data, -1)
        return data.reshape(self.mesh.shape + data.shape[1:])
        
    def unshape_cells(self, data):
        shape = data.shape
        if shape[-len(self.mesh.shape):] == self.mesh.shape:
            shape = shape[:-len(self.mesh.shape)] + (-1,)
        return data.reshape(shape)
        
    @property
    def mesh(self):
        if not hasattr(self, "_mesh") or self._mesh is None:
            topology = Topology.from_xml(document=self.document, node=self.node.getElementsByTagName("Topology")[0])
            geometry = Geometry.from_xml(document=self.document, node=self.node.getElementsByTagName("Geometry")[0])
            
            nx, ny, nz = topology.nxnynz
            dx, dy, dz = geometry.dxdydz
            from fipy.meshes.uniformGrid3D import UniformGrid3D
            self._mesh = UniformGrid3D(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz) + geometry.origin
            
        return self._mesh
    
class Topology(_Node):
    @classmethod
    def from_xml(cls, document, node):
        topologyType = node.getAttribute("TopologyType") or node.getAttribute("Type")
        if topologyType == "3DCoRectMesh":
            return ThreeDCoRectMeshTopology(document=document, node=node)
        else:
            raise Exception("Unknown TopologyType: ", topologyType)
        
class StructuredTopology(Topology):
    @property
    def nxnynz(self):
        return xml_to_array(self.node.getAttribute("Dimensions"))[::-1] - 1
        
class ThreeDCoRectMeshTopology(StructuredTopology):
    pass
    
class Geometry(_Node):
    @classmethod
    def from_xml(cls, document, node):
        geometryType = node.getAttribute("GeometryType") or node.getAttribute("Type")
        if geometryType == "ORIGIN_DXDYDZ":
            return ORIGIN_DXDYDZGeometry(document=document, node=node)
        else:
            raise Exception("Unknown GeometryType: ", geometryType)

class ORIGIN_DXDYDZGeometry(Geometry):
    @property
    def origin(self):
        o = DataItem.from_xml(document=self.document, 
                              node=self.node.getElementsByTagName("DataItem")[0]).array
        return o[numerix.newaxis, ...]
        
    @property
    def dxdydz(self):
        return DataItem.from_xml(document=self.document, 
                                 node=self.node.getElementsByTagName("DataItem")[1]).array
        
class CollectionGrid(Grid):
    @classmethod
    def from_name(cls, document, name=""):
        collection = document.node.createElement("Grid")
        collection.setAttribute("Name", name)
        collection.setAttribute("GridType", "Collection")
        
        return cls(document=document, node=collection)
        
    @property
    def variables(self):
        vars = []
        for node in self.node.childNodes:
            if node.nodeType == minidom.Node.ELEMENT_NODE:
                if node.tagName == "Grid":
                    vars.extend(Grid.from_xml(document=self.document, node=node).variables)
                elif node.tagName == "Attribute":
                    vars.append(Attribute(document=self.document, node=node).variable)
        
        return vars
