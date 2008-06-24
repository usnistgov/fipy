__docformat__ = 'restructuredtext'

from mayaviViewer import MayaviViewer
from mayaviSurfactantViewer import MayaviSurfactantViewer

__all__ = ["MayaviViewer", "MayaviSurfactantViewer"]

def make(*args, **kwargs):
    import warnings
    warnings.warn("'MayaviViewer' should be used instead of 'make'", DeprecationWarning, stacklevel=2)
    return MayaviViewer(*args, **kwargs)
