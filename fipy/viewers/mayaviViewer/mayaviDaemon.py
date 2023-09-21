"""A simple script that polls a data file for changes and then updates
the Mayavi pipeline automatically.

This script is based heavily on the `poll_file.py` example in the Mayavi distribution.


This script is to be run like so::

 $ mayavi2 -x mayaviDaemon.py <options>

Or::

 $ python mayaviDaemon.py <options>

Run::

 $ python mayaviDaemon.py --help

to see available options.
"""
from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

# Standard imports.
import os
import signal
import sys

# Enthought library imports
try:
    from mayavi.plugins.app import Mayavi
    from mayavi.sources.vtk_file_reader import VTKFileReader
    from pyface.timer.api import Timer
    from mayavi import mlab
    from tvtk.api import tvtk
except ImportError as e:
    from enthought.mayavi.plugins.app import Mayavi
    from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
    from enthought.pyface.timer.api import Timer
    from enthought.mayavi import mlab

# FiPy library imports
from fipy.tools.numerix import array, concatenate, where, zeros

__all__ = ["MayaviDaemon", "main"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

######################################################################
class MayaviDaemon(Mayavi):
    """Given a file name and a mayavi2 data reader object, this class
    polls the file for any changes and automatically updates the
    mayavi pipeline.
    """

    _viewers = []

    def parse_command_line(self, argv):
        """Parse command line options.

        Parameters
        ----------
        argv : :obj:`list` of :obj:`str`
            The command line arguments
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

        parser.add_option("--fps", action="store", dest="fps", type="float", default=1.0,
                          help="frames per second to attempt to display")

        (options, args) = parser.parse_args(argv)

        self._lockfname = options.lock
        self._cellfname = options.cell
        self._facefname = options.face
        self._bounds = [options.xmin, options.xmax,
                        options.ymin, options.ymax,
                        options.zmin, options.zmax]

        self._datamin = options.datamin
        self._datamax = options.datamax

        self._fps = options.fps

    @staticmethod
    def _examine_data(source, datatype, bounds):
        """Determine contents of source

        Parameters
        ----------
        source : tvtk.DataSet
        datatype : str
            either "cell_data" or "point_data"
        bounds : array_like
            boundaries of existing data sets

        Returns
        -------
        has : dict
            whether each rank is present in data set
        bounds : array_like
            boundaries of data sets
        """
        ranks = ["scalars", "vectors", "tensors"]
        has = dict((rank, False) for rank in ranks)

        if source is not None:
            # Newer versions of mayavi (> 4.7?) store AssignAttribute objects
            # in outputs, so this clumsy bit is to extract the underlying
            # DataSet objects.
            # This is a clear sign that we're using this completely wrong,
            # but, eh, who cares?
            sourceoutputs = [out if isinstance(out, tvtk.DataSet)
                             else out.trait_get()['output']
                             for out in source.outputs]

            for rank in ranks:
                tmp = [out.trait_get()[datatype].trait_get()[rank]
                       for out in sourceoutputs]
                tmp = [out for out in tmp if out is not None]
                has[rank] = (len(tmp) > 0)

            bounds = concatenate((bounds,
                                  [out.bounds for out in sourceoutputs]),
                                 axis=0)

        return has, bounds

    def run(self):
        MayaviDaemon._viewers.append(self)

        mlab.clf()

        self.cellsource = self.setup_source(self._cellfname)
        self.has_cell, bounds = self._examine_data(source=self.cellsource, datatype="cell_data",
                                                   bounds=zeros((0, 6), 'l'))

        self.facesource = self.setup_source(self._facefname)
        self.has_face, bounds = self._examine_data(source=self.facesource, datatype="point_data",
                                                   bounds=bounds)

        boundsmin = bounds.min(axis=0)
        boundsmax = bounds.max(axis=0)

        bounds = (boundsmin[0], boundsmax[1],
                  boundsmin[2], boundsmax[3],
                  boundsmin[4], boundsmax[5])

        self._bounds = where(self._bounds == array((None,)),
                             bounds,
                             self._bounds).astype(float)

        self.view_data()

        # Poll the lock file.
        self.timer = Timer(1000 / self._fps, self.poll_file)

    def __del__(self):
        dir = None
        for fname in [self._cellfname, self._facefname, self._lockfname]:
            if fname and os.path.isfile(fname):
                os.unlink(fname)
                if not dir:
                    dir = os.path.dirname(fname)
        if dir:
            os.rmdir(dir)

    @staticmethod
    def _sigint_handler(signum, frame):
        for viewer in MayaviDaemon._viewers:
            viewer.__del__()
        raise SystemExit("MayaviDaemon cleaned up")

    def poll_file(self):
        if os.path.isfile(self._lockfname):
            self.update_pipeline(self.cellsource)
            self.update_pipeline(self.facesource)
            with open(self._lockfname, 'r') as lock:
                filename = lock.read()
            if len(filename) > 0:
                mlab.savefig(filename)
            os.unlink(self._lockfname)

    def update_pipeline(self, source):
        """Override this to do something else if needed.
        """
        if source is not None:
            source.scene.disable_render = True
            source.scene.anti_aliasing_frames = 0
            # Force the reader to re-read the file.
            source.reader.modified()
            source.update()
            # Propagate the changes in the pipeline.
            source.data_changed = True
            source.scene.disable_render = False

    def setup_source(self, fname):
        """Given a VTK file name `fname`, this creates a mayavi2 reader
        for it and adds it to the pipeline.  It returns the reader
        created.
        """
        if fname is None:
            return None

        source = VTKFileReader()
        source.initialize(fname)
        mlab.pipeline.add_dataset(source)

        return source

    def clip_data(self, src):
        if hasattr(mlab.pipeline, "data_set_clipper"):
            clip = mlab.pipeline.data_set_clipper(src)
            clip.filter.inside_out = True

            clip.widget.widget_mode = 'Box'
            clip.widget.widget.place_factor = 1.
            clip.widget.widget.place_widget(self._bounds)
            clip.widget.update_implicit_function()

            clip.widget.visible = False
        else:
            import warnings
            warnings.warn("Mayavi r24017 or newer needed for data_set_clipper()", UserWarning, stacklevel=2)
            clip = src

        return clip

    def _view_data(self, source, has, has_scale_bar, cell_data=False):
        """Determine contents of source

        Parameters
        ----------
        source : tvtk.DataSet
        has : dict
            whether each rank is present in data set
        has_scale_bar : bool
            whether a scale bar has already been created
        cell_data : bool
            whether source contains cell_data that may need conversion
            to point_data

        Returns
        -------
        has_scale_bar : bool
            whether a scale bar has been created
        """
        if source is not None:
            clip = self.clip_data(source)

            if has["scalars"]:
                s = mlab.pipeline.surface(clip, vmin=self._datamin, vmax=self._datamax)
                if not has_scale_bar:
                    s.module_manager.scalar_lut_manager.show_scalar_bar = True
                    has_scale_bar = True
            if cell_data:
                clip = mlab.pipeline.cell_to_point_data(clip)
            if has["vectors"]:
                v = mlab.pipeline.vectors(clip, vmin=self._datamin, vmax=self._datamax)
                if not has_scale_bar:
                    v.module_manager.scalar_lut_manager.show_scalar_bar = True
                    has_scale_bar = True

        return has_scale_bar

    def view_data(self):
        """Sets up the mayavi pipeline for the visualization.
        """
        has_scale_bar = self._view_data(source=self.cellsource, has=self.has_cell,
                                        has_scale_bar=False, cell_data=True)
        has_scale_bar = self._view_data(source=self.facesource, has=self.has_face,
                                        has_scale_bar=has_scale_bar)

signal.signal(signal.SIGINT, MayaviDaemon._sigint_handler)
try:
    signal.signal(signal.SIGHUP, MayaviDaemon._sigint_handler)
except AttributeError:
    # not available on Windows
    pass
signal.signal(signal.SIGTERM, MayaviDaemon._sigint_handler)

def main(argv=None):
    """Simple helper to start up the mayavi application.

    This returns the running application.
    """
    m = MayaviDaemon()
    m.main(argv)
    return m

if __name__ == '__main__':
    main(sys.argv[1:])
