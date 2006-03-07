#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mayaviViewer.py"
 #                                    created: 9/14/04 {2:48:25 PM} 
 #                                last update: 3/3/06 {4:07:06 PM} { 2:45:36 PM}
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##


__docformat__ = 'restructuredtext'

import mayavi
import pyvtk

from fipy.viewers.viewer import Viewer

class MayaviViewer(Viewer):
    """
    The `MayaviViewer` creates viewers with the Mayavi_ python
    plotting package.

    .. _Mayavi: http://mayavi.sourceforge.net/

    Issues with the `MayaviViewer` are

      - `_getOrderedCellVertexIDs()` doesn't return the correct ordering
        for 3D meshes.  This may be okay for tets and wedges but will
        break for hexahedrons.

      - Different element types can not be displayed for 3D
        meshes. This is an ordering issue for the cell data. Could get
        round this either by implementing a method such as
        `var.getVertexVariable()` and use point data, or reordering the
        variable data via [tets, wedges, hexs] and keep using cell
        data. First option is cleaner. Second option is less work.

      - Should this class be split into various dimensions? Is it
        useful to display data with different dimension is same viewer?

    """
        
    def __init__(self, vars, limits = None, title = None):
        """
        Create a `MayaviViewer`.
        
        :Parameters:

          - `vars`: a `CellVariable` or tuple of `CellVariable` objects to plot
          - `limits`: a dictionary with possible keys `xmin`, `xmax`,
            `ymin`, `ymax`, `zmin`, `zmax`, `datamin`, `datamax`.  A 1D
            Viewer will only use `xmin` and `xmax`, a 2D viewer will also
            use `ymin` and `ymax`, and so on.  All viewers will use
            `datamin` and `datamax`.  Any limit set to a (default) value of
            `None` will autoscale.
          - `title`: displayed at the top of the Viewer window

        """

        Viewer.__init__(self, vars = vars, limits = limits, title = title)

        self._viewer = mayavi.mayavi()

        self.structures = []
            
        for var in self.vars:
            self.structures.append(self._getStructure(var.getMesh()))
                                       
    def _getStructure(self, mesh):

        cellVertexIDs = mesh._getOrderedCellVertexIDs()

        import fipy.tools.numerix as numerix
        lengths = len(cellVertexIDs[0]) - numerix.sum(numerix.MA.getmaskarray(cellVertexIDs), index = 1)
        
        cellDict = {2 : [], 4: [], 6: [], 8: [], 'polygon' : []}

        cellVertexListIDs = [list(cellVertexIDs[i][:lengths[i]]) for i in range(mesh.getNumberOfCells())]

        if mesh.getDim() == 2:
            cellDict['polygon'] = cellVertexListIDs
        else:

            if mesh.getDim() == 3:
                print "Warning: The Matayvi viewer may or may not render 3D meshes correctly"
                print "The method Mesh._getOrderedCellVertexIDs() needs to be fixed to return the"
                print "correct ordering for 3D meshes that fits with pyvtk paradigm."

            import sets
            if len(sets.Set(lengths)) > 1:
                raise TypeError, 'The MayaviViewer can only render 3D data with cells of the same type'
            else:
                if lengths[0] in cellDict.keys():
                    cellDict[lengths[0]] = cellVertexListIDs
                else:
                    raise TypeError, 'These cell types can not be rendered by the MayaviViewer'


                
        coords = numerix.zeros((mesh.getVertexCoords().shape[0], 3), 'd')
        coords[:,:mesh.getDim()] = mesh.getVertexCoords()

        return pyvtk.UnstructuredGrid(points = coords,
                                      line = cellDict[2],
                                      tetra = cellDict[4],
                                      wedge = cellDict[6],
                                      voxel = cellDict[8],
                                      polygon = cellDict['polygon'])

    def plot(self, filename = None):
        """
        Plot the `CellVariable` as a contour plot.

        :Parameters:
          - `filename`: The name of the file for PNG hard copies.
        
        """

        import os
        import tempfile
        
        for var, structure in zip(self.vars, self.structures):
            name = var.getName()
            if name is '':
                name = 'default'

            celldata = pyvtk.CellData(pyvtk.Scalars(var[:], name = name, lookup_table = 'default'))
            data = pyvtk.VtkData(structure, "mydata", celldata)
            
            (f, fileName) = tempfile.mkstemp('.vtk')
            data.tofile(fileName)
            self._viewer.open_vtk(fileName, config=0)
            
            os.close(f)
            os.remove(fileName)
            
            self._viewer.load_module('SurfaceMap', 0)
            rw = self._viewer.get_render_window()
            rw.z_plus_view()

            ## display legend
            dvm = self._viewer.get_current_dvm()
            mm = dvm.get_current_module_mgr()
            slh = mm.get_scalar_lut_handler()
            slh.legend_on.set(1)
            slh.legend_on_off()

            ## display legen with correct range
            slh.range_on_var.set(1)
            slh.v_range_on_var.set(1)

            xmax = self._getLimit('datamax')
            if xmax is None:
                xmax = numerix.max(var)

            xmin = self._getLimit('datamin')
            if xmin is None:
                xmin = numerix.min(var)
            
            slh.range_var.set((xmin, xmax))
            slh.set_range_var()

            slh.v_range_var.set((numerix.min(var), numerix.max(var)))
            slh.set_v_range_var()

            self._viewer.Render()

        if filename is not None:
            self._viewer.renwin.save_png(filename)


if __name__ == '__main__':
    from fipy.meshes.grid1D import Grid1D
    from fipy.variables.cellVariable import CellVariable
    vars = [CellVariable(value = range(3), mesh = Grid1D(nx = 3, dx = 1.))]

    from fipy.meshes.tri2D import Tri2D
    triMesh = Tri2D()
    from fipy.meshes.grid2D import Grid2D
    gridMesh = Grid2D(nx = 3)
    gridMesh += (1, 0)

    compositeMesh = gridMesh + triMesh
    compositeMesh += (0, 2)
    vars += [CellVariable(value = range(7), mesh = compositeMesh)]
    
    from fipy.meshes.grid3D import Grid3D
    mesh3D = Grid3D(nx = 3)
    mesh3D += (0, 0, 2)
    vars += [CellVariable(value = range(3), mesh = mesh3D)]

    viewer = MayaviViewer(vars)
    viewer.plot()
    raw_input("finished")
    
    
