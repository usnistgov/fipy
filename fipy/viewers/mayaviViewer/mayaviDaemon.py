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
from mmap import mmap
import os
from os.path import join, abspath, dirname
import sys

# Enthought library imports
from enthought.mayavi.scripts import mayavi2
from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
from enthought.mayavi.modules.outline import Outline
from enthought.mayavi.modules.contour_grid_plane import ContourGridPlane
from enthought.mayavi.filters.cell_to_point_data import CellToPointData
from enthought.pyface.timer.api import Timer


######################################################################
# `Pollster` class.
class Pollster(object):
    """Given a file name and a mayavi2 data reader object, this class
    polls the file for any changes and automatically updates the
    mayavi pipeline.
    """
    def __init__(self, fname, data):
        """Initialize the object.

        Parameters:
        -----------
        fname -- filename to poll.
        data -- the MayaVi source object to update.
        """
        self.fname = fname
        f = os.open(self.fname, os.O_RDWR)
        self.mmap = mmap(f, 1024)
        self.mmap.write("READY")
        self.mmap.seek(0)

        self.data = data
#         self.last_stat = os.stat(fname)
        
#     def poll_file(self):
#         # Check the file's time stamp.
#         s = os.stat(self.fname)
#         if s.st_mtime == self.last_stat.st_mtime:
#             return
#         else:
#             self.last_stat = s
#             self.update_pipeline()

    def poll_file(self):
        msg = self.mmap.readline()
        self.mmap.seek(0)
        if msg.startswith("PLOT"):
            self.update_pipeline()
            self.mmap.write("READY")
            self.mmap.seek(0)

    def update_pipeline(self):
        """Override this to do something else if needed.
        """
        # Force the reader to re-read the file.
        d = self.data
        d.reader.modified()
        d.update()
        # Propagate the changes in the pipeline.
        d.data_changed = True


                       
def setup_data(fname):
    """Given a VTK file name `fname`, this creates a mayavi2 reader
    for it and adds it to the pipeline.  It returns the reader
    created.
    """
    # 'mayavi' is always defined on the interpreter.
    mayavi.new_scene()
    d = VTKFileReader()
    d.initialize(fname)
    s = mayavi.add_source(d)
    return s, d

def view_data(src):
    """Sets up the mayavi pipeline for the visualization.
    """
    # 'mayavi' is always defined on the interpreter.
    o = Outline()
    mayavi.add_module(o)
    
#     c = ContourGridPlane()
#     mayavi.add_module(c)
#     c.grid_plane.position = 16
#     c.module_manager.scalar_lut_manager.show_scalar_bar = True

    from enthought.mayavi import mlab
    s = mlab.pipeline.surface(src) #,extent=extent,vmin=datamin,vmax=datamax)
#     s.scene.anti_aliasing_frames = 0
    s.module_manager.scalar_lut_manager.show_scalar_bar = True

@mayavi2.standalone
def main():
    # Change this to suit your needs.  Edit the file after running this
    # script and the pipeline should be updated automatically.

    vtkfname = sys.argv[1]
    mmapfname = sys.argv[2]
    
    
    src, data = setup_data(vtkfname)
    view_data(src)

    # Poll the file.
    p = Pollster(mmapfname, data)
    timer = Timer(1000, p.poll_file)

    mayavi.fipytimer = timer
    
    # To stop polling the file do:
    #timer.Stop()

if __name__ == '__main__':
    main()

