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

# Based on poll_file.py
#
# Author: Prabhu Ramachandran <prabhu@aero.iitb.ac.in>
# Copyright (c) 2006-2007, Enthought Inc.
# License: BSD Style.

# Standard imports.
import os
import sys

# Enthought library imports
from enthought.mayavi.plugins.app import Mayavi
from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
from enthought.pyface.timer.api import Timer
from enthought.mayavi import mlab

# FiPy library imports
from fipy.tools.numerix import array, concatenate, where, zeros

######################################################################
class MayaviDaemon(Mayavi):
    """Given a file name and a mayavi2 data reader object, this class
    polls the file for any changes and automatically updates the
    mayavi pipeline.
    """
    def parse_command_line(self, argv):
        """Parse command line options.

        Parameters
        ----------

        - argv : `list` of `strings`

          The list of command line arguments.
        """
        from optparse import OptionParser
        usage = "usage: %prog [options]"
        parser = OptionParser(usage)

        parser.add_option("-l", "--lock", action="store", dest="lock", type="string", default=None,
                          help="path of lock file")

        parser.add_option("-c", "--cell", action="store", dest="cell", type="string", default=None,
                          help="path of cell vtk file")

        parser.add_option("-f", "--face", action="store", dest="face", type="string", default=None,
                          help="path of face vtk file")

        parser.add_option("--xmin", action="store", dest="xmin", type="float", default=None,
                          help="minimum x value")

        parser.add_option("--xmax", action="store", dest="xmax", type="float", default=None,
                          help="maximum x value")

        parser.add_option("--ymin", action="store", dest="ymin", type="float", default=None,
                          help="minimum y value")

        parser.add_option("--ymax", action="store", dest="ymax", type="float", default=None,
                          help="maximum y value")

        parser.add_option("--zmin", action="store", dest="zmin", type="float", default=None,
                          help="minimum z value")

        parser.add_option("--zmax", action="store", dest="zmax", type="float", default=None,
                          help="maximum z value")

        parser.add_option("--datamin", action="store", dest="datamin", type="float", default=None,
                          help="minimum data value")

        parser.add_option("--datamax", action="store", dest="datamax", type="float", default=None,
                          help="maximum data value")

        (options, args) = parser.parse_args(argv)
        
        self.lockfname = options.lock
        self.cellfname = options.cell
        self.facefname = options.face
        self.bounds = [options.xmin, options.xmax, 
                       options.ymin, options.ymax, 
                       options.zmin, options.zmax]
                       
        self.datamin = options.datamin
        self.datamax = options.datamax
        
    def run(self):
        # 'mayavi' is always defined on the interpreter.
        mayavi.new_scene()

        bounds = zeros((0, 6))
        
        self.cellsource = self.setup_source(self.cellfname)
        if self.cellsource is not None:
            tmp = [out.cell_data.scalars for out in self.cellsource.outputs \
                   if out.cell_data.scalars is not None]
            self.has_cell_scalars = (len(tmp) > 0)
            tmp = [out.cell_data.vectors for out in self.cellsource.outputs \
                   if out.cell_data.vectors is not None]
            self.has_cell_vectors = (len(tmp) > 0)
            tmp = [out.cell_data.tensors for out in self.cellsource.outputs \
                   if out.cell_data.tensors is not None]
            self.has_cell_tensors = (len(tmp) > 0)

            bounds = concatenate((bounds, 
                                  [out.bounds for out in self.cellsource.outputs]),
                                 axis=0)


        self.facesource = self.setup_source(self.facefname)
        if self.facesource is not None:
            tmp = [out.point_data.scalars for out in self.facesource.outputs \
                   if out.point_data.scalars is not None]
            self.has_face_scalars = (len(tmp) > 0)
            tmp = [out.point_data.vectors for out in self.facesource.outputs \
                   if out.point_data.vectors is not None]
            self.has_face_vectors = (len(tmp) > 0)
            tmp = [out.point_data.tensors for out in self.facesource.outputs \
                   if out.point_data.tensors is not None]
            self.has_face_tensors = (len(tmp) > 0)
            
            bounds = concatenate((bounds, 
                                  [out.bounds for out in self.facesource.outputs]),
                                 axis=0)
                                 
        boundsmin = bounds.min(axis=0)
        boundsmax = bounds.max(axis=0)
        
        bounds = (boundsmin[0], boundsmax[1], 
                  boundsmin[2], boundsmax[3], 
                  boundsmin[4], boundsmax[5])

        self.bounds = where(self.bounds == array((None,)),
                            bounds, 
                            self.bounds).astype(float)

        self.view_data()

        # Poll the lock file.
        self.timer = Timer(1000, self.poll_file)

    def poll_file(self):
        if os.path.isfile(self.lockfname):
            self.update_pipeline(self.cellsource)
            self.update_pipeline(self.facesource)
            lock = file(self.lockfname, 'r')
            filename = lock.read()
            lock.close()
            if len(filename) > 0:
                mlab.savefig(filename)
            os.unlink(self.lockfname)

    def update_pipeline(self, source):
        """Override this to do something else if needed.
        """
        if source is not None:
            # Force the reader to re-read the file.
            source.reader.modified()
            source.update()
            # Propagate the changes in the pipeline.
            source.data_changed = True
        
    def setup_source(self, fname):
        """Given a VTK file name `fname`, this creates a mayavi2 reader
        for it and adds it to the pipeline.  It returns the reader
        created.
        """
        if fname is None:
            return None
            
        source = VTKFileReader()
        source.initialize(fname)
        mayavi.add_source(source)
        
        return source
        
    def clip_data(self, src):
        clip = mlab.pipeline.data_set_clipper(self.cellsource)
        clip.filter.inside_out = True

        clip.widget.widget_mode = 'Box'
        clip.widget.widget.place_factor = 1.
        clip.widget.widget.place_widget(self.bounds)
        clip.widget.update_implicit_function()

        clip.widget.visible = False
        
        return clip

    def view_data(self):
        """Sets up the mayavi pipeline for the visualization.
        """
        from enthought.tvtk.api import tvtk
            
        has_scale_bar = False
        if self.cellsource is not None:
            clip = self.clip_data(self.cellsource)
            
            if self.has_cell_scalars:
                s = mlab.pipeline.surface(clip, vmin=self.datamin, vmax=self.datamax)
                s.module_manager.scalar_lut_manager.show_scalar_bar = True
                has_scale_bar = True
            p = mlab.pipeline.cell_to_point_data(clip)
            if self.has_cell_vectors:
                v = mlab.pipeline.vectors(p, vmin=self.datamin, vmax=self.datamax)
                if not has_scale_bar:
                    v.module_manager.scalar_lut_manager.show_scalar_bar = True
                    has_scale_bar = True

        if self.facesource is not None:
            clip = self.clip_data(self.facesource)

            if self.has_face_scalars:
                s = mlab.pipeline.surface(clip, vmin=self.datamin, vmax=self.datamax)
                if not has_scale_bar:
                    s.module_manager.scalar_lut_manager.show_scalar_bar = True
                    has_scale_bar = True
            if self.has_face_vectors:
                v = mlab.pipeline.vectors(clip, vmin=self.datamin, vmax=self.datamax)
                if not has_scale_bar:
                    v.module_manager.scalar_lut_manager.show_scalar_bar = True
                    has_scale_bar = True

def main(argv=None):
    """Simple helper to start up the mayavi application.  This returns
    the running application."""
    m = MayaviDaemon()
    m.main(argv)
    return m

if __name__ == '__main__':
    main(sys.argv[1:])
