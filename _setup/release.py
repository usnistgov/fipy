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

    def _build_unix_distribution(self):
        """Create Unix source distribution"""

        self._remove_manifest()
        shutil.copyfile("MANIFEST-UNIX.in", "MANIFEST.in")
        run_setup("setup.py", ["sdist"])
        os.remove("MANIFEST.in")

    def _build_windows_distribution(self):
        """Create Windows source distribution

        Contains executable installer and examples"""

        import versioneer

        version = versioneer.get_version()

        self._remove_manifest()

        shutil.copyfile("MANIFEST-WINDOWS.in", "MANIFEST.in")
        run_setup("setup.py", ["sdist", "--dist-dir=dist-windows", "--formats=zip"])
        shutil.move(
            os.path.join("dist-windows", "FiPy-{}.zip".format(version)),
            os.path.join("dist", "FiPy-{}.win32.zip".format(version)),
        )
        os.rmdir("dist-windows")
        os.remove("MANIFEST.in")

    def run(self):
        if self.unix:
            self._build_unix_distribution()
        if self.windows:
            self._build_windows_distribution()
