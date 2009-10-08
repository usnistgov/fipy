#!/usr/bin/env python
"""A simple script that polls a data file for changes and then updates
the mayavi pipeline automatically.

This script is based heavily on the poll_file.py exampe in the mayavi distribution.


This script is to be run like so::

 $ mayavi2 -x mayaviDaemon.py ???

Or::

 $ python mayaviDaemon.py ???
 
The script currently defaults to using the example data in
examples/data/heart.vtk.  You can try editing that data file or change
this script to point to other data which you can edit.
"""

# Author: Jonathan Guyer <guyer@nist.gov>

# Author: Prabhu Ramachandran <prabhu@aero.iitb.ac.in>
# Copyright (c) 2006-2007, Enthought Inc.
# License: BSD Style.

# Standard imports.
import os
import sys

# Enthought library imports
from enthought.mayavi.scripts import mayavi2
from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
from enthought.pyface.timer.api import Timer


######################################################################
class MayaviDaemon(object):
    """Given a file name and a mayavi2 data reader object, this class
    polls the file for any changes and automatically updates the
    mayavi pipeline.
    """
    def __init__(self, cellfname, facefname, lockfname):
        """Initialize the object.

        Parameters:
        -----------
        cellfname -- VTK file holding cell data.
        facefname -- VTK file holding face data.
        lockfname -- filename to poll.
        """
        self.cellfname = cellfname
        self.facefname = facefname
        self.lockfname = lockfname
        
        # 'mayavi' is always defined on the interpreter.
        mayavi.new_scene()

        self.cellsource = self.setup_source(self.cellfname)
        self.facesource = self.setup_source(self.facefname)
        self.view_data()

        # Poll the lock file.
        self.timer = Timer(1000, self.poll_file)


    def poll_file(self):
        if os.path.isfile(self.lockfname):
            self.update_pipeline(self.cellsource)
            self.update_pipeline(self.facesource)
            os.unlink(self.lockfname)

    def update_pipeline(self, data):
        """Override this to do something else if needed.
        """
        # Force the reader to re-read the file.
        data.reader.modified()
        data.update()
        # Propagate the changes in the pipeline.
        data.data_changed = True
        
    def setup_source(self, fname):
        """Given a VTK file name `fname`, this creates a mayavi2 reader
        for it and adds it to the pipeline.  It returns the reader
        created.
        """
        source = VTKFileReader()
        source.initialize(fname)
        mayavi.add_source(source)
        
        return source

    def view_data(self):
        """Sets up the mayavi pipeline for the visualization.
        """
        from enthought.mayavi import mlab
            
        o = mlab.pipeline.outline(self.cellsource)
        scalars = [out.cell_data.scalars for out in self.cellsource.outputs \
                   if out.cell_data.scalars is not None]
        if len(scalars) > 0:
            s = mlab.pipeline.surface(self.cellsource) #,extent=extent,vmin=datamin,vmax=datamax)
#             s.module_manager.scalar_lut_manager.show_scalar_bar = True
        p = mlab.pipeline.cell_to_point_data(self.cellsource)
        vectors = [out.cell_data.vectors for out in self.cellsource.outputs \
                   if out.cell_data.vectors is not None]
        if len(vectors) > 0:
            v = mlab.pipeline.vectors(p)

        scalars = [out.point_data.scalars for o in self.facesource.outputs \
                   if out.point_data.scalars is not None]
        if len(scalars) > 0:
            s = mlab.pipeline.surface(self.facesource) #,extent=extent,vmin=datamin,vmax=datamax)
    #     s.module_manager.scalar_lut_manager.show_scalar_bar = True
        p = mlab.pipeline.cell_to_point_data(self.facesource)
        vectors = [out.point_data.vectors for out in self.facesource.outputs \
                   if out.point_data.vectors is not None]
        if len(vectors) > 0:
            v = mlab.pipeline.vectors(p)

@mayavi2.standalone
def main():
    # Change this to suit your needs.  Edit the file after running this
    # script and the pipeline should be updated automatically.

    cellfname = sys.argv[1]
    facefname = sys.argv[2]
    lockfname = sys.argv[3]
    
    mayavi.fipydaemon = MayaviDaemon(cellfname, facefname, lockfname)
    
    # To stop polling the file do:
    #mayavi.fipydaemon.timer.Stop()

if __name__ == '__main__':
    main()

