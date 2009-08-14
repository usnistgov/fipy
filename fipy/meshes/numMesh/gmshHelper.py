class GmshHelper(object):

	def __init__(self,version):

                #Names describing the columns in the 'Elements' section of the gmsh file
		colNames=[]
		
		colNames.append({
		'ElemNumber':0,
		'ElemType':1,
		'PhysicalNum':2,
		'PartitionNum':-1,
		'VertexStart':5
		})

		colNames.append({
		'ElemNumber':0,
		'ElemType':1,
		'PhysicalNum':2,
		'PartitionNum':5,
		'VertexStart':6
		})

		colNames=colNames[version-1]

                #The types of primitives used by gmsh. Because we don't use the subdivisions, each categeroy contains many different primitives of the same shape.
		types= {
		'point':[15],
		'line':[1,8,26,27,28],
		'triangle':[2,9,20,21,22,23,24,25],
		'quadrangle':[3,10,16],
		'tetrahedron':[4,11,29,30,31],
		'hexahedron':[5,12,17],
		'triangular prism':[6,13,18],
		'square pyramid':[7,14,19],
		}

                #The names of each type of primitive.
		invTypes = {}
		for k,v in types.iteritems():
			for num in v:
				invTypes[num]=k

                #The number of dimensions each type of primitive lies in, by name
		_dims = {
		'point':0,
		'line':1,
		'triangle':2,
		'quadrangle':2,
		'tetrahedron':3,
		'hexahedron':3,
		'triangular prism':3,
		'square pyramid':3
		}

                #The number of dimensions each type of primitive lies in, by number
		dims = {}
                #The list of primitives, by number, that lie in each dimension
		invDims = {}
		for k,v in _dims.iteritems():
			for n in types[k]:
				dims[n]=v
			if v in invDims.keys():
				invDims[v].extend(types[k])
				invDims[v].sort()
			else:
				invDims[v]=types[k]

                #The number of vertices each type of primitive contains
		numVerts = {1:2, 2:3, 3:4, 4:4, 5:8, 6:6, 7:5, 8:3, 9:6, 10:9, 11:10, 12:27, 13:18, 14:14, 15:1, 16:8, 17:20, 18:15, 19:13, 20:9, 21:10, 22:12, 23:15, 24:15, 25:21, 26:4, 27:5, 28:6, 29:20, 30:35, 31:56}

		from numpy import ma,zeros,array

                #A masked array describing the relationship between the cells and faces for each primitive, by name
		faces = {
		'point':ma.MaskedArray(array([]),array([]),dtype='int').reshape((0,0)),
		'line':ma.MaskedArray(array([0,1]),zeros(2),dtype='int').reshape((2,1)),
		'triangle':ma.MaskedArray(array([0,1,1,2,2,0]),zeros(6),dtype='int').reshape((3,2)),
		'quadrangle':ma.MaskedArray(array([0,1,1,2,2,3,3,0]),zeros(8),dtype='int').reshape((4,2)),
		'tetrahedron':ma.MaskedArray(array([0,1,2,0,1,3,0,2,3,1,2,3]),zeros(12),dtype='int').reshape((4,3)),
		'hexahedron':ma.MaskedArray(array([0,1,2,3,0,1,5,4,0,3,7,4,1,2,6,5,2,3,7,6,4,5,6,7]),zeros(24),dtype='int').reshape((6,4)),
		'triangular prism':ma.MaskedArray(array([0,1,2,0,0,1,4,3,0,2,5,3,1,2,5,4,3,4,5,0]),array([0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]),dtype='int').reshape((5,4)),
		'square pyramid':ma.MaskedArray(array([0,1,2,3,0,1,4,0,0,3,4,0,1,2,4,0,2,3,4,0]),array([0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1]),dtype='int').reshape((5,4))
		}

		for arr in faces.values():
			arr.fill_value=0

		self.colNames=colNames
		self.types=types
		self.invTypes=invTypes
		self.dims=dims
		self.invDims=invDims
		self.numVerts=numVerts
		self.faces=faces

	#Gets the array associated with a number rather than a name
	def getArr(self,t):
		return self.faces[self.invTypes[t]]
