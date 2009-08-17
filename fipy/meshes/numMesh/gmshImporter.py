#!/usr/bin/env python

##
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "gmshImport.py"
 #
 #  Author: Alexander Mont <alexander.mont@nist.gov>
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
 #

__docformat__ = 'restructuredtext'

r"""

This module takes a Gmsh output file (`.msh`) and converts it into a
FiPy mesh. 

Gmsh generates unstructured meshes, which may contain a significant
amount of non-orthogonality and it is very difficult to directly
control the amount of non-orthogonality simply by manipulating Gmsh
parameters. Therefore, it is necessary to take into account the
possibility of errors arising due to the non-orthogonality of the
mesh. To test the degree of error, an experiment was conducted using a
simple 1D diffusion problem with constant diffusion coefficients and
boundary conditions as follows: fixed value of 0 on the left side,
fixed value of 1 on the right side, and a fixed flux of 0 on the top
and bottom sides. The analytical solution is clearly a uniform
gradient going from left to right. this problem was implemented using
a Cartesian Grid2D mesh with each interior vertex displaced a short
distance in a random direction to create non-orthogonality. Then the
root-mean-square error was plotted against the root-mean-square
non-orthogonality. The error in each cell was calculated by simply
subtracting the analytical solution at each cell center from the
calculated value for that cell. The non-orthogonality in each cell is
the average, weighted by face area, of the sines of the angles between
the face normals and the line segments joining the cells. Thus, the
non-orthogonality of a cell can range from 0 (every face is orthogonal
to its corresponding cell-to-cell line segment) to 1 (only possible in
a degenerate case). This test was run using 500 separate 20x20 meshes
and 500 separate 10x10 meshes, each with the interior vertices moved
different amounts so as to created different levels of
non-orthogonality. The results are shown below.

Results for 20x20 mesh:

.. image:: fipy/meshes/numMesh/orthoerrorgraph.pdf
   :height: 100px
   :width: 200px

Results for 10x10 mesh:

.. image:: fipy/meshes/numMesh/orthoerrorcoarsegraph.pdf
    :height: 100px
    :width: 200px

It is clear from the graphs that finer meshes decrease the error due
to non-orthogonality, and that even with a reasonably coarse mesh the
error is quite low. However, note that this test is only for a simple
1D diffusion problem with a constant diffusion coefficient, and it is
unknown whether the results will be significantly different with more
complicated problems.

Test cases:

   >>> newmesh = GmshImporter3D('fipy/meshes/numMesh/testgmsh.msh')
   >>> print newmesh.getVertexCoords()
   [[ 0.   0.5  1.   0.5  0.5]
    [ 0.   0.5  0.   1.   0.5]
    [ 0.   1.   0.   0.   0.5]]

   >>> print newmesh._getFaceVertexIDs()
   [[0 0 0 0 0 0 1 1 1 2]
    [1 1 1 2 2 3 2 2 3 3]
    [2 3 4 3 4 4 3 4 4 4]]

   >>> print newmesh._getCellFaceIDs()
   [[0 1 3 6]
    [2 2 4 7]
    [4 5 5 8]
    [7 8 9 9]]

   >>> mesh = GmshImporter2DIn3DSpace('fipy/meshes/numMesh/GmshTest2D.msh')
   >>> print mesh.getVertexCoords()
   [[ 0.   1.   0.5  0.   1.   0.5  0.   1. ]
    [ 0.   0.   0.5  1.   1.   1.5  2.   2. ]
    [ 0.   0.   0.   0.   0.   0.   0.   0. ]]

   >>> mesh = GmshImporter2D('fipy/meshes/numMesh/GmshTest2D.msh')
   >>> print mesh.getVertexCoords()
   [[ 0.   1.   0.5  0.   1.   0.5  0.   1. ]
    [ 0.   0.   0.5  1.   1.   1.5  2.   2. ]]

   >>> print mesh._getFaceVertexIDs()
   [[ 1 0 3 2 4 2 4 3 5 6 4 4 5 7 6]
    [ 0 2 0 1 1 3 2 4 3 3 5 7 6 5 7]]
   
   >>> print mesh._getCellFaceIDs()
   [[ 1  1  3  7  7  8 13 14]
    [ 3  5  6  6 10 12 10 13]
    [ 0  2  4  5  8  9 11 12]]
   
The following test case is to test the handedness of the mesh to check
it does not return negative volumes. Firstly we set up a list with
tuples of strings to be read by gmsh. The list provide instuctions to
gmsh to form a circular mesh.

   >>> cellSize = 0.7
   >>> radius = 1.
   >>> lines = ['cellSize = ' + str(cellSize) + ';\n',
   ...           'radius = ' + str(radius) + ';\n',
   ...           'Point(1) = {0, 0, 0, cellSize};\n',
   ...           'Point(2) = {-radius, 0, 0, cellSize};\n',
   ...           'Point(3) = {0, radius, 0, cellSize};\n',
   ...           'Point(4) = {radius, 0, 0, cellSize};\n',
   ...           'Point(5) = {0, -radius, 0, cellSize};\n',
   ...           'Circle(6) = {2, 1, 3};\n',
   ...           'Circle(7) = {3, 1, 4};\n',
   ...           'Circle(8) = {4, 1, 5};\n',
   ...           'Circle(9) = {5, 1, 2};\n',
   ...           'Line Loop(10) = {6, 7, 8, 9} ;\n',
   ...           'Plane Surface(11) = {10};\n']
   >>> lines = ''.join(lines)

Check that the sign of the mesh volumes is correct

   >>> mesh = GmshImporter2D(lines)
   >>> print mesh.getCellVolumes()[0] > 0
   1
      
Reverse the handedness of the mesh and check the sign

   >>> lines[7:12] = ['Circle(6) = {3, 1, 2};\n',
   ...                'Circle(7) = {4, 1, 3};\n',
   ...                'Circle(8) = {5, 1, 4};\n',
   ...                'Circle(9) = {2, 1, 5};\n',
   ...                'Line Loop(10) = {9, 8, 7, 6};\n',]

   >>> mesh = GmshImporter2D(lines)
   >>> print mesh.getCellVolumes()[0] > 0
   1
   
"""

from gmshHelper import GmshHelper

def genMshFile(filename,dimensions):
	'''Generates a .msh file using gmsh.  Can be passed in a piece of code as a string, a .gmsh filename, a .geo filenae, or a .msh filename.
	   If the filename ends in .msh, simply returns the same filename'''
	import tempfile,os
	rmGmshFile=False
	ext = filename[filename.rindex('.'):]
	if ext == '.msh':
		return filename
	elif ext not in ['.gmsh','.geo']:
		(f,name) = tempfile.mkstemp('.gmsh','fipy',os.getcwd(),text=True)
		os.write(f,filename)
		os.close(f)
		filename=name
		rmGmshFile=True
	(f,name2) = tempfile.mkstemp('.msh','fipy',os.getcwd())
	os.system('gmsh -'+str(dimensions)+' -v 0 -format msh -o '+name2+' '+filename)
	if rmGmshFile:
		os.remove(name)
	os.close(f)
	return name2
		

def formatMshFile(filename,newName):
	'''This method formats a file with name `filename such that it can be used by numpy.genfromtxt.
	   A file formatted by this function is required to use the "genMeshFromMshFile" method'''
	import string
	newName = newName or filename[:-4]+'CSV.msh'
	#Create a temp file in which all the whitespace is replaced by commas
	#Also find the maximum number of columns in each section of the file, for use later
	fi = open(filename,'rU')
	tmp = open("tmp.msh",'w')
	tags = {} #Dictionary from a tag ($Elements) to the numbermaximum number of columns between that tag and the previous one.  This means that in a properly structured gmsh file, any tag not starting with $End will be associated with 0.
	maxCols = 0
	currLine = fi.readline()
	prevTag=''
	while currLine:
		if (currLine[0]=='$'):
			tag=currLine.strip()[1:]
			tags[tag] = maxCols
			maxCols = 0
			if prevTag:
				tags[tag]=(tags[tag],prevTag)
				tags[prevTag]=(tags[prevTag],tag)
				prevTag=''
			else:
				prevTag=tag
				
		currLine = string.join(currLine.split(),',')
		maxCols = max(maxCols,currLine.count(',')+1)
		tmp.write(currLine+'\n')
		currLine = fi.readline()
	fi.close()
	tmp.flush()
	tmp.close()

	#Open the new file to write to, and make it so that each line between two tags has the same amount of commas.
	fi = open("tmp.msh",'rU')
	tmp = open(newName,'w')
	
	currLine = fi.readline()
	inTag = '' #This will be the name of the block we are currently in.  If we just hit an end tag, it will be zero.
	currCols=0
	while currLine:
		numCommas=0
		if (currLine[0]=='$'):
			if inTag:
				inTag = ''
			else:
				#The maximum number of columns is stored in the $End tag corresponding to this one.
				inTag = currLine.strip()[1:]
				currCols = tags[tags[inTag][1]][0]
		else:
			numCommas = currCols-1-currLine.count(',')
		tmp.write(currLine.strip())
		if numCommas > 0:
	       		tmp.write(','*numCommas)
		tmp.write('\n')
		currLine = fi.readline()
	
	fi.close()
	tmp.flush()
	tmp.close()
	
	import os
	os.remove('tmp.msh')
	
	return newName

def findTags(filename):
	'''This method returns a dictionary from tags in a .msh file to the position of the character the tag starts at and the position it ends at.
	   It adds on the matching tag, IE it matches the first tag with the second and vice versa, the third with the fourth and so on.
	   If the file is correctly structured, this should match tags with end tags.
	   For example, if the first line were "$Elements", the dictionary would make "Elements" to (0,10) since the \n character is included in this count.
	   This is useful for parsing the files, as you can just read from the end of one tag to the beginning of its corresponding End tag.'''
	tags = {}
	fi = open(filename,'rU')
	prevPos = fi.tell()
	currLine = fi.readline()
	currPos = fi.tell()
	prevTag = ''
	while currLine:
		if currLine[0]=='$':
			tag = currLine.strip()[1:]
			tags[tag]=(prevPos,currPos)
			if prevTag:
				tags[prevTag]=tags[prevTag]+(tag,)
				tags[tag]=tags[tag]+(prevTag,)
				prevTag=''
			else:
				prevTag=tag
		prevPos = currPos
		currLine = fi.readline()
		currPos = fi.tell()
	fi.close()
	return tags
	

def subFile(fi,start,end):
	'''This method returns a file-like object that reads only between the two given points in the file.'''
	from StringIO import StringIO
	fi.seek(start)
	return StringIO(fi.read(end-start))
	
def invIDs(ids):
	from numpy.ma import zeros,ones,arange,indices,MaskedArray,bitwise_or
	if type(ids) != type(MaskedArray(0)):
		ids = MaskedArray(ids,zeros(ids.shape),dtype=int)
	l = ids.max()+1
	u = arange(l)
	u.mask = zeros(u.shape)
	e = bitwise_or.reduce(u == ids[...,None])
	ind = indices(ids.shape)[1,0][...,None]
	m = e.sum(axis=0).max()
	z = MaskedArray(zeros((m,l)),ones((m,l)),dtype=int)

	for n in xrange(len(u)):
		opp = (e*ind)[e[:,n],n]
		z.mask[:opp.size,n]=False
		z[:opp.size,n]=opp
	return z
	
def GmshImporter(filename,dimensions=3,shapeDim=None,formatted=False,keepFile=False,csvFile=None):
	'''This method generates a FiPy mesh using gmsh.
	   :Parameters:
	       filename: the name of the file to be parsed (.msh,.gmsh, or .geo) or a string of commands in the gmsh language
	       dimensions: the number of dimensions in which the mesh is embedded.
	         This will simply use the first n coordinates of the points.
	         Default: 3
	       shapeDim: the maximum dimensionality of the shapes to be generated
	         (a triangle can lie in 3 space but still only has 2 dimensions)
		 Default: dimensions
	       formatted: whether or the the ".msh" file is in the format requred by this method.
	         Having the file pref-formatted will save a few seconds.
		 To pre-format it, just call the "formatMshFile" method from another python script
		   and use the file that ends with "...CSV.msh" as the filename.
	         Default: False
	       keepFile: whether or not to keep the .msh file generated by this method (if any)
	         Default: False
	       csvFile: The filename to output the CSV file to.  This is only useful if keepFile is True.
	         Default: None'''
	shapeDim = shapeDim or dimensions
	
	import os
	
	#format the file, get the new filename
	MSHfilename = genMshFile(filename,shapeDim)
	if MSHfilename != filename:
		formatted = False
	if not formatted:
		CSVfilename = formatMshFile(MSHfilename,csvFile)
	else:
		CSVfilename=MSHfilename
	if MSHfilename != filename:
		os.remove(MSHfilename)
	#Find the tags in the .msh file
	tags = findTags(CSVfilename)

	if tags.has_key('MeshFormat'):
		version = 2
	else:
		version = 1

	elms = (version - 1) and 'Elements' or 'ELM'
	endElms = (version - 1) and 'EndElements' or 'ENDELM'
	nodes = (version - 1) and 'Nodes' or 'NOD'
	endNodes = (version - 1) and 'EndNodes' or 'ENDNOD'

	gmshHelper = GmshHelper(version)

	from numpy import genfromtxt,ma,unique1d,arange,zeros,ones,sort,array,indices

	#Load the lines under the Elements tag into an array (Ignore the first line as it is just the number of elements)
	fi = open(CSVfilename)
	
	arr = genfromtxt(subFile(fi,tags[elms][1],tags[endElms][0]),delimiter=',',usemask=True,dtype='int')[1:]
	types=unique1d(arr[:,gmshHelper.colNames['ElemType']]) #The types of primitives found in the .msh file
	maxDims=0
	maxFacesPerCell=0
	maxVerticesPerFace=0
	for t in types: #Loop throught the types to find the maximum dimension, number of faces per cell, and number of vertices per face
		maxDims=max(maxDims,gmshHelper.dims[t])
		maxFacesPerCell=max(maxFacesPerCell,gmshHelper.getArr(t).shape[0])
		maxVerticesPerFace=max(maxVerticesPerFace,gmshHelper.getArr(t).shape[1])
	if dimensions < maxDims: #Trying to load tetrahedrons into 2D space... doesn't really work.  However, loading triangles in 3D space does, so we let the other inequality pass through.
		raise ValueError('You cannot have %d dimensional cells in a %d dimensional mesh'%(maxDims,dimensions))
	cellFaceIDs=ma.MaskedArray(zeros((0,maxFacesPerCell)),zeros((0,maxFacesPerCell)),dtype='int')
	faceVertexIDs=ma.MaskedArray(zeros((0,maxVerticesPerFace)),zeros((0,maxVerticesPerFace)),dtype='int')
	totFaces=0
	maxDimTypes=[]
	for t in types: #loop through the types of cells of the highest dimension in the .msh file
		if gmshHelper.dims[t]==maxDims:
			maxDimTypes.append(t)
			currArr=gmshHelper.getArr(t)
			#Get all the elements of this type
			elems=arr[arr[:,gmshHelper.colNames['ElemType']]==t,gmshHelper.colNames['VertexStart']:gmshHelper.colNames['VertexStart']+gmshHelper.numVerts[t]]
			#Construct Cell to Face IDs for all the cells of this type
			cellFacePart=ma.MaskedArray(arange(totFaces,totFaces+(currArr.shape[0]*elems.shape[0])),zeros(currArr.shape[0]*elems.shape[0]),dtype='int').reshape((elems.shape[0],currArr.shape[0]))
			#Use the array in gmshHelper to give us the FaceVertex IDs for all the cells of this type
			faceVertexPart=elems[:,currArr]
			faceVertexPart.mask[:]|=currArr.mask
			faceVertexPart = faceVertexPart.reshape(-1,currArr.shape[1])
			#Fill in masked values if there are more Faces per Cell in some other type of cell
			cellFacePart=ma.concatenate((cellFacePart,ma.MaskedArray(zeros((cellFacePart.shape[0],maxFacesPerCell-cellFacePart.shape[1])),ones((cellFacePart.shape[0],maxFacesPerCell-cellFacePart.shape[1])),dtype='int')),axis=1)
			#Fill in masked values if there are more Vertices per Face in some other type of cell
			faceVertexPart=ma.concatenate((faceVertexPart,ma.MaskedArray(zeros((faceVertexPart.shape[0],maxVerticesPerFace-faceVertexPart.shape[1])),ones((faceVertexPart.shape[0],maxVerticesPerFace-faceVertexPart.shape[1])),dtype='int')),axis=1)
			#Add the information about cells of this type to the overall array
			cellFaceIDs=ma.concatenate((cellFaceIDs,cellFacePart),axis=0)
			faceVertexIDs=ma.concatenate((faceVertexIDs,faceVertexPart),axis=0)

	#concatenation sometimes meses up the masks if they are all false, so lets fix that
	if cellFaceIDs.mask.shape != cellFaceIDs.shape:
		tmp = zeros(cellFaceIDs.shape)
		tmp[:]=cellFaceIDs.mask
		cellFaceIDs.mask=tmp
	if faceVertexIDs.mask.shape != faceVertexIDs.shape:
		tmp = zeros(faceVertexIDs.shape)
		tmp[:]=faceVertexIDs.mask
		faceVertexIDs.mask=tmp

	#Create a view of the Face to Vertex IDs that is 1D, with the vertices
	#    of each face sorted into ascending order.
	#This allows us to use numpy.unique1d to find out what faces are duplicated, and where.
	faceVertexIDs.fillValue=-1
	cellFaceIDs.fillValue=-1
	sorted,ind,inv=unique1d(sort(faceVertexIDs.filled(),axis=1).view('S%d'%(faceVertexIDs.itemsize*faceVertexIDs.size/faceVertexIDs.shape[0])),return_index=True,return_inverse=True)
	#Update all the Cell to Face IDs to the new positions of the faces
	cellFaceIDs[cellFaceIDs.mask==False]=inv[cellFaceIDs[cellFaceIDs.mask==False]]
	cellFaceIDs=cellFaceIDs.swapaxes(0,1)
	#Update the Face to Vertex IDs with the new order
	faceVertexIDs=faceVertexIDs[ind].swapaxes(0,1)
	#Get the Cell to Vertex IDs (since this is currently unused in fipy, no reason to spend the time doing it.)
	#cellVertexIDs=arr[(arr[None,:,gmshHelper.colNames['ElemType']]==array(maxDimTypes)[:,None]).sum(axis=0).astype('bool'),gmshHelper.colNames['VertexStart']:].swapaxes(0,1)

	#Get the vertex coordinates
	vertexCoords = genfromtxt(subFile(fi,tags[nodes][1],tags[endNodes][0]),delimiter=',')[1:,:dimensions+1].swapaxes(0,1)
	
	vertexPointIDs=zeros((2,vertexCoords.shape[1]),dtype='int')
	vertexPointIDs[0]=indices((vertexCoords.shape[1],))[0]
	vertexPointIDs[1]=vertexCoords[0]
	vertexPointIDs=ma.MaskedArray(vertexPointIDs,zeros(vertexPointIDs.shape),dtype='int')
	vertexCoords=vertexCoords[1:]

	faceVertexIDs=ma.MaskedArray(((faceVertexIDs==vertexPointIDs[1,:,None,None]).data*vertexPointIDs[0,:,None,None]).sum(axis=0),faceVertexIDs.mask)
	fi.close()
	if not keepFile:
		os.remove(CSVfilename)
	
	faceVertexIDs.fillValue=-1
	cellFaceIDs.fillValue=-1
	
	#Order the vertices in the FaceVertexIDs 

	#create the mesh
	from mesh import Mesh
	from mesh1D import Mesh1D
	from mesh2D import Mesh2D
	if maxDims==1:
		cla = Mesh1D
	elif maxDims==2:
		cla = Mesh2D
	else:
		cla = Mesh
	return cla(vertexCoords,faceVertexIDs.filled(),cellFaceIDs.filled())

def GmshImporter2D(filename,formatted=False,keepFile=False,csvFile=None):
	return GmshImporter(filename,2,2,formatted,keepFile,csvFile)

def GmshImporter3D(filename,formatted=False,keepFile=False,csvFile=None):
	return GmshImporter(filename,3,3,formatted,keepFile,csvFile)

def GmshImporter2DIn3DSpace(filename,formatted=False,keepFile=False,csvFile=None):
	return GmshImporter(filename,3,2,formatted,keepFile,csvFile)
