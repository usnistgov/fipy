"""setuptools command to prepare FiPy for release"""

from distutils.core import Command
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
    user_options = [("version=", None, "Version string to tag and assign to release")]

    def initialize_options(self):
        import versioneer

        self.version = versioneer.get_version()

    def finalize_options(self):
        pass

    def _remove_manifest(self):
        """Remove MANIFEST file

        probably no longer needed, MANIFEST was ancient history?"""

        try:
            os.remove("MANIFEST")
        except OSError, _:
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

        self._remove_manifest()
        run_setup("setup.py", ["bdist", "--formats=wininst"])

        self._remove_manifest()
        fname = "FiPy-{}.win32.exe".format(self.version)
        os.symlink(os.path.join("dist", fname), fname)
        shutil.copyfile("MANIFEST-WINDOWS.in", "MANIFEST.in")
        run_setup("setup.py", ["sdist", "--dist-dir=dist-windows", "--formats=zip"])
        os.unlink(fname)
        shutil.move(
            os.path.join("dist-windows", "FiPy-{}.zip".format(self.version)),
            os.path.join("dist", "FiPy-{}.win32.zip".format(self.version)),
        )
        os.rmdir("dist-windows")
        os.remove("MANIFEST.in")

    def run(self):
        run_setup("setup.py", ["bdist_egg"])
        run_setup("setup.py", ["build_docs", "--pdf", "--html", "--cathartic"])

        self._build_unix_distribution()
        self._build_windows_distribution()
