#!/usr/bin/env python
## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "time.py"
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

from fipy.io.xdmf.grid import CollectionGrid, Grid
from fipy.io.xdmf.node import _Node

class TimeSeries(_Node):
    @classmethod
    def empty(cls, document):
        series = document.node.createElement("Grid")
        series.setAttribute("Name", "TimeSeries")
        series.setAttribute("GridType", "Collection")
        series.setAttribute("CollectionType", "Temporal")
        
        return cls(document=document, node=series)

    def __setitem__(self, index, value):
        # find and remove existing index?
        collection = CollectionGrid.from_name(document=self.document, name="Step " + index)
        collection += Time.from_value(document=self.document, value=index)
        if not isinstance(value, (list, tuple)):
            value = (value,)
            
        h5filename = "%s-%d.h5" % (os.path.splitext(self.document.filename)[0], len(self.node.childNodes))

        for node in NodesFromValues(document=self.document, values=value, h5filename=h5filename):
            collection += node
            
        if os.path.exists(h5filename):
            # if a corresponding .h5 file has been created, time stamp it,
            # otherise the time stored in the .xmf is enough
            import h5py
            f = h5py.File(h5filename, mode='a')
            f.parent.attrs["Time"] = index
            f.close()
        
#         self.node.replaceChild
#         self += GridFromValues(doc=self.doc, values=value)
        self += collection
        
        self.document.save()
        
    @property
    def _nodes(self):
        return [node for node in self.node.childNodes 
                if node.nodeType == minidom.Node.ELEMENT_NODE 
                and node.tagName == "Grid"
                and len(node.getElementsByTagName("Time")) > 0]
                     
    def __getitem__(self, index):
        nodes = [node for node in self._nodes 
                 if (node.getElementsByTagName("Time")[0].getAttribute("Value") == index)]
        
        grid = Grid.from_xml(document=self.document, node=nodes[0])
        
        return grid.variables
        
    def keys(self):
        return [node.getElementsByTagName("Time")[0].getAttribute("Value") for node in self._nodes]


class Time(_Node):
    @classmethod
    def from_value(cls, document, value):
        time = document.node.createElement("Time")
        time.setAttribute("Value", value)

        return cls(document=document, node=time)

def NodesFromValues(document, values, h5filename):
    from fipy.variables.meshVariable import _MeshVariable

    nodes = []
    for value in values:
        if isinstance(value, _MeshVariable):
            if len(nodes) > 0 and isinstance(nodes[-1], Grid) and nodes[-1].mesh is value.mesh:
                grid = nodes[-1]
            else:
                grid = value.mesh._to_xdmf(document=document)
                nodes.append(grid)
            
            grid += value._to_xdmf(document=document, grid=grid, h5filename=h5filename)
        else:
            nodes.append(value._to_xdmf(document=document, grid=None, h5filename=h5filename))
            
    return nodes
    
