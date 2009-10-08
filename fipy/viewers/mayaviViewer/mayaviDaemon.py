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
from enthought.pyface.timer.api import Timer


######################################################################
# `Pollster` class.
class Pollster(object):
    """Given a file name and a mayavi2 data reader object, this class
    polls the file for any changes and automatically updates the
    mayavi pipeline.
    """
    def __init__(self, fname, cd, fd):
        """Initialize the object.

        Parameters:
        -----------
        fname -- filename to poll.
        cd -- the MayaVi cell source object to update.
        fd -- the MayaVi face source object to update.
        """
        self.fname = fname
        self.celldata = cd
        self.facedata = fd

    def poll_file(self):
        if os.path.isfile(self.fname):
            self.update_pipeline()
            os.unlink(self.fname)

    def update_pipeline(self):
        """Override this to do something else if needed.
        """
        # Force the reader to re-read the file.
        cd = self.celldata
        cd.reader.modified()
        cd.update()
        # Propagate the changes in the pipeline.
        cd.data_changed = True

        # Force the reader to re-read the file.
        fd = self.facedata
        fd.reader.modified()
        fd.update()
        # Propagate the changes in the pipeline.
        fd.data_changed = True


                       
def setup_data(cellfname, facefname):
    """Given a VTK file name `fname`, this creates a mayavi2 reader
    for it and adds it to the pipeline.  It returns the reader
    created.
    """
    # 'mayavi' is always defined on the interpreter.
    mayavi.new_scene()
    
    cd = VTKFileReader()
    cd.initialize(cellfname)
    cs = mayavi.add_source(cd)
    
    fd = VTKFileReader()
    fd.initialize(facefname)
    fs = mayavi.add_source(fd)
    
    return cs, cd, fs, fd

def view_data(cs, fs):
    """Sets up the mayavi pipeline for the visualization.
    """
    # 'mayavi' is always defined on the interpreter.
#     o = Outline()
#     mayavi.add_module(o)
    
    from enthought.mayavi import mlab
    o = mlab.pipeline.outline(cs)
    s = mlab.pipeline.surface(cs) #,extent=extent,vmin=datamin,vmax=datamax)
    s.module_manager.scalar_lut_manager.show_scalar_bar = True
    p = mlab.pipeline.cell_to_point_data(cs)
    v = mlab.pipeline.vectors(p)

    s = mlab.pipeline.surface(fs) #,extent=extent,vmin=datamin,vmax=datamax)
#     s.module_manager.scalar_lut_manager.show_scalar_bar = True
    p = mlab.pipeline.cell_to_point_data(fs)
    v = mlab.pipeline.vectors(p)

@mayavi2.standalone
def main():
    # Change this to suit your needs.  Edit the file after running this
    # script and the pipeline should be updated automatically.

    vtkcellfname = sys.argv[1]
    vtkfacefname = sys.argv[2]
    vtklockfname = sys.argv[3]
    
    
    cs, cd, fs, fd = setup_data(vtkcellfname, vtkfacefname)
    view_data(cs, fs)

    # Poll the file.
    p = Pollster(vtklockfname, cd, fd)
    timer = Timer(1000, p.poll_file)

    mayavi.fipytimer = timer
    
    # To stop polling the file do:
    #timer.Stop()

if __name__ == '__main__':
    main()

