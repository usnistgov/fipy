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

from fipy.io.xdmf.node import _Node

class _Document(_Node):
    def __init__(self, filename, mode, doc):
        _Node.__init__(self, document=self, node=doc)
        
        self.filename = filename
        self.mode = mode
        
        xdmf = doc.getElementsByTagName("Xdmf")[0]
        self.domain = xdmf.getElementsByTagName("Domain")[0]
        
    def save(self):
        f = open(self.filename, mode=self.mode)
        self.node.writexml(f, indent="    ", addindent="    ", newl="\n")
        f.close()
        
def Open(filename, mode='r'):
    if mode.startswith('w'):
        imp = minidom.getDOMImplementation('')
        dt = imp.createDocumentType("Xdmf", None, "Xdmf.dtd")
        doc = imp.createDocument(None, "Xdmf", dt)
        xdmf = doc.documentElement
        xdmf.setAttribute("xmlns:xi", "http://www.w3.org/2003/XInclude")
        
        domain = doc.createElement("Domain")
        xdmf.appendChild(domain)
    else:
        doc = minidom.parse(filename)

    return _Document(filename=filename, mode=mode, doc=doc)
        
def GridFromValues(doc, values):
    """
    >= 1 _MeshVariable with same Mesh
    >= 1 _MeshVariable with different Mesh
    >= 1 Variable not _MeshVariable, plus >= 1 _MeshVariable with same Mesh
    >= 1 Variable not _MeshVariable, plus >= 1 _MeshVariable with different Mesh
    
    tuple coming out should be in same order as tuple going in, even if means redeclaring meshes
    """
    from fipy.variables.meshVariable import _MeshVariable
    
    unmeshed = []
    meshed = []
    meshes = []
    for value in values:
        if isinstance(value, _MeshVariable):
            meshed.append(value)
            if value.mesh not in meshes:
                meshes.append(value.mesh)
        else:
            unmeshed.append(value)
            
    if len(meshes) == 0:
        pass
    elif len(meshes) == 1:
        pass
    else:
        grid = CollectionGrid(doc)
#         for value in values:
#             grid += 
            
#     meshless = [value for value in values if not isinstance(value, _MeshVariable)]
#     values = [value for value in values if isinstance(value, _MeshVariable)]
#     meshed = {}
#     while len(values) > 0:
#         mesh = values[0].mesh
#         meshed[id(mesh)] = [value for value in values if value.mesh is mesh]
#         values = [value for value in values if value.mesh is not mesh]
        
    return meshes
    
if __name__ == "__main__":
    from fipy.io.xdmf import Open, TimeSeries
    
    xdmf = Open("test.xmf", "w")
    ts = TimeSeries.empty(document=xdmf)
    xdmf.domain.appendChild(ts.node)
    
    from fipy import *
    mesh1 = Grid3D(nx=2, ny=3, nz=4)
    x, y, z = mesh1.cellCenters
    var1 = Variable(name="var1", value=3.0)
    var2 = CellVariable(mesh=mesh1, name="var2", value=x*y*z)
    varx = CellVariable(mesh=mesh1, name="x", value=x)
    vary = CellVariable(mesh=mesh1, name="y", value=y)
    varz = CellVariable(mesh=mesh1, name="z", value=z)

    mesh2 = Grid3D(nx=3, ny=2, nz=14, dz=0.5) + [[3], [1], [1]]
    x, y, z = mesh2.cellCenters
    var3 = CellVariable(mesh=mesh2, value=x*y*z)

#     meshes = GridFromValues(doc=xdmf.doc, values=(var1, var2, var3))
    
#     print "meshes:", meshes
    
#     print "meshless:", meshless
#     print "meshed:", meshed
    
    for time in arange(0, 1, 0.1):
        var2.value = var2 + time
        ts["%0.1g" % time] = (var1, var2, varx, vary, varz) #, var3)
        xdmf.save()
        
    xdmf2 = Open("test.xmf", "r")
    ts = TimeSeries(document=xdmf2, node=xdmf2.domain.getElementsByTagName("Grid")[0])
#     var1a, var2a, varxa, varya, varza, var2grad = ts["0.9"]
    var1a, var2a, varxa, varya, varza = ts["0.9"]
    
    print var1a.allclose(var1)
    print var2a.allclose(var2)
    print varxa.allclose(varx)
    print varya.allclose(vary)
    print varza.allclose(varz)
#     print var2grad.allclose(var2.grad)

#     print xdmf
#     xdmf.save()
