import os
import platform
import subprocess
import sys
from xml.dom.minidom import Document

import fipy

__all__ = ["Vitals"]

class Vitals(Document):
    """Returns XML formatted information about current FiPy environment
    """

    def __init__(self):
        Document.__init__(self)

        self.top = self.createElementNS("http://www.ctcms.nist.gov/fipy", "FiPy")
        Document.appendChild(self, self.top)

        self.appendChild(self.tupleToXML(sys.argv, "sys.argv"))

        version = self.createElement("version")
        self.appendChild(version)
        version.appendChild(self.createTextNode(fipy.__version__))

        path = self.createElement("path")
        fipypath = os.path.dirname(fipy.__file__)
        path.appendChild(self.createTextNode(fipypath))
        self.appendChild(path)

        self.appendChild(self.svn(fipypath))

        self.appendChild(self.dictToXML(os.environ, "environ"))

        pyth = self.createElement("python")

        implementation = self.createElement("implementation")
        if hasattr(platform, "python_implementation"):
            implementation.appendChild(self.createTextNode(platform.python_implementation()))
        else:
            implementation.appendChild(self.createTextNode("unknown"))
        pyth.appendChild(implementation)

        pversion = self.createElement("version")
        pversion.appendChild(self.createTextNode(platform.python_version()))
        pyth.appendChild(pversion)

        compiler = self.createElement("compiler")
        compiler.appendChild(self.createTextNode(platform.python_compiler()))
        pyth.appendChild(compiler)

        pyth.appendChild(self.tupleToXML(platform.python_build(), "build",
                                         keys=("buildno", "builddate")))

        pyth.appendChild(self.tupleToXML(platform.architecture(), "architecture",
                                         keys=("bits", "linkage")))

        pyth.appendChild(self.tupleToXML(platform.uname(), "uname",
                                         keys=("system", "node", "release", "version", "machine", "processor")))

        self.appendChild(pyth)

    def appendChild(self, child):
        self.top.appendChild(child)

    def dictToXML(self, d, name):
        elem = self.createElement(name)
        for key, value in d.items():
            keyelem = self.createElement(key)
            keyelem.appendChild(self.createTextNode(str(value)))
            elem.appendChild(keyelem)

        return elem

    def tupleToXML(self, t, name, keys=None):
        elem = self.createElement(name)
        if keys is not None:
            for key, value in zip(keys, t):
                keyelem = self.createElement(key)
                keyelem.appendChild(self.createTextNode(str(value)))
                elem.appendChild(keyelem)
        else:
            for value in t:
                elem.appendChild(self.createTextNode(str(value)))

        return elem

    def svncmd(self, cmd, *args):
        elem = self.createElement(cmd)
        p = subprocess.Popen(["svn", cmd] + list(args), stdout=subprocess.PIPE)
        info = p.communicate()[0]
        if p.returncode == 0:
            elem.appendChild(self.createTextNode(info))
        else:
            raise OSError

        return elem

    def svn(self, *args):
        elem = self.createElement("svn")
        for cmd in ["info", "status", "diff"]:
            try:
                elem.appendChild(self.svncmd(cmd, *args))
            except:
                pass

        return elem

    def __str__(self):
        return self.toprettyxml()

    def save(self, fname):
        f = open(fname, 'w')
        self.writexml(f, indent="    ", addindent="    ", newl="\n")
        f.close()

    def appendInfo(self, name, svnpath=None, **kwargs):
        """append some additional information, possibly about a project under a separate svn repository
        """
        elem = self.dictToXML(kwargs, name)
        if svnpath is not None:
            elem.appendChild(self.svn(svnpath))
        self.appendChild(elem)

if __name__ == "__main__":
    v = Vitals()

    solar = v.createElement("solar")
    solar.appendChild(v.svn("/Users/guyer/Documents/research/codes/solar-dimensionless"))
    v.appendChild(solar)

    print v
