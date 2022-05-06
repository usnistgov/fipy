"""setuptools command to prepare FiPy for release"""

from distutils.core import Command
import glob
import os
import shutil

from setuptools.sandbox import run_setup

__all__ = ["release"]


class release(Command):
    """Prepare FiPy for release

    Generates tarball and excutable Windows installer"""

    description = "Prepare the FiPy release artifacts"

    # List of option tuples: long name, short name (None if no short
    # name), and help string.
    user_options = [('unix', None, "create a tarball source distribution"),
                    ('windows', None, "create an executable installer for MS Windows"),
                    ('all', None, "create unix and Windows distributions"),
                   ]

    def initialize_options(self):
        self.unix = 0
        self.windows = 0
        self.all = 0

    def finalize_options(self):
        if self.all:
            self.unix = 1
            self.windows = 1

    def _remove_manifest(self):
        """Remove MANIFEST file

        probably no longer needed, MANIFEST was ancient history?"""

        try:
            os.remove("MANIFEST")
        except OSError as _:
            pass

    def _build_source_distribution(self, formats):
        """Create source distribution"""

        self._remove_manifest()
        run_setup("setup.py", ["sdist", "--formats=" + ",".join(formats)])

    def run(self):
        formats = []
        if self.unix:
            formats.append("gztar")
        if self.windows:
            formats.append("zip")
        self._build_source_distribution(formats)
