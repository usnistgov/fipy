__docformat__ = 'restructuredtext'

from mayaviViewer import MayaviViewer

__all__ = ["MayaviViewer"]

def make(*args, **kwargs):
    import warnings
    warnings.warn("'MayaviViewer' should be used instead of 'make'", DeprecationWarning, stacklevel=2)
    return MayaviViewer(*args, **kwargs)
