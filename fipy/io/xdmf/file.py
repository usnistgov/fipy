import os
from xml.dom import minidom

from fipy.tools import numerix

class _Node(object):
    def __init__(self, document, node):
        self.document = document
        self.node = node
        
    def __str__(self):
        return self.node.toprettyxml()
        
    def __iadd__(self, value):
        self.node.appendChild(value.node)
        return self
        
    def __getitem__(self, index):
        return self.node.getElementsByTagName(index)
        

class Document(_Node):
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

    return Document(filename=filename, mode=mode, doc=doc)
        
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
        collection = CollectionGrid.from_name(document=self.document, name="Step %g" % index)
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
        
    def __getitem__(self, index):
        nodes = [node for node in self.node.childNodes 
                 if node.nodeType == minidom.Node.ELEMENT_NODE 
                 and node.tagName == "Grid"]
        node = [node for node in nodes 
                if (node.getElementsByTagName("Time")[0].getAttribute("Value") == str(index))][0]
        
        grid = Grid.from_xml(document=self.document, node=node)
        
#         raise Exception("stop!")

        return grid.variables

class Time(_Node):
    @classmethod
    def from_value(cls, document, value):
        time = document.node.createElement("Time")
        time.setAttribute("Value", str(value))

        return cls(document=document, node=time)

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

def _xml_to_array(xml):
    from fipy.tools import numerix
    from StringIO import StringIO
    return numerix.loadtxt(StringIO(xml)).ravel()

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
        
def NodesFromValues(document, values, h5filename):
    from fipy.variables.meshVariable import _MeshVariable
    from fipy.variables import CellVariable, FaceVariable

    nodes = []
    for value in values:
        if isinstance(value, _MeshVariable):
            if len(nodes) > 0 and isinstance(nodes[-1], Grid) and nodes[-1].mesh is value.mesh:
                grid = nodes[-1]
            else:
                grid = MeshGrid.from_Mesh(document=document, mesh=value.mesh)
                nodes.append(grid)
            
            if isinstance(value, CellVariable):
                attr = CellAttribute.from_CellVariable(document=document, var=value, grid=grid, h5filename=h5filename)
            elif isinstance(value, FaceVariable):
                attr = FaceAttribute.from_FaceVariable(document=document, var=value, grid=grid, h5filename=h5filename)

            grid += attr
        else:
            nodes.append(Attribute.from_Variable(document=document, var=value))
            
    return nodes
    
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
        
    @classmethod
    def from_Mesh(cls, document, mesh):
        from fipy.meshes.uniformGrid3D import UniformGrid3D
        if isinstance(mesh, UniformGrid3D):
            grid = ThreeDCoRectMeshGrid.from_Mesh(document=document, mesh=mesh)
        else:
            raise Exception("Don't know how!")
            
        return grid
        
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
# the following are properly specific Mesh methods?
    @classmethod
    def from_Mesh(cls, document, mesh):
        grid = document.node.createElement("Grid")
        grid.setAttribute("Name", mesh.__class__.__name__)
        grid.setAttribute("GridType", "Uniform")
        
        topology = document.node.createElement("Topology")
        topology.setAttribute("TopologyType", "3DCoRectMesh")
        shape = list(mesh.shape)
        shape.reverse()
        topology.setAttribute("Dimensions", " ".join(str(v+1) for v in shape))
        grid.appendChild(topology)
        
        geometry = document.node.createElement("Geometry")
        geometry.setAttribute("GeometryType", "ORIGIN_DXDYDZ")
        grid.appendChild(geometry)

        origin = XMLDataItem.from_array(document=document, 
                                        name="Origin",
                                        arr=mesh.origin.ravel())
        geometry.appendChild(origin.node)
        
        spacing = XMLDataItem.from_array(document=document, 
                                         name="Spacing",
                                         arr=(mesh.dx, mesh.dy, mesh.dz))
        geometry.appendChild(spacing.node)

        return cls(document=document, node=grid, mesh=mesh)
        
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
        return _xml_to_array(self.node.getAttribute("Dimensions"))[::-1] - 1
        
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



if __name__ == "__main__":
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
        ts[time] = (var1, var2, varx, vary, varz) #, var3)
        xdmf.save()
        
    xdmf2 = Open("test.xmf", "r")
    ts = TimeSeries(document=xdmf, node=xdmf2.domain.getElementsByTagName("Grid")[0])
#     var1a, var2a, varxa, varya, varza, var2grad = ts[0.9]
    var1a, var2a, varxa, varya, varza = ts[0.9]
    
    print var1a.allclose(var1)
    print var2a.allclose(var2)
    print varxa.allclose(varx)
    print varya.allclose(vary)
    print varza.allclose(varz)
    print var2grad.allclose(var2.grad)

#     print xdmf
#     xdmf.save()
