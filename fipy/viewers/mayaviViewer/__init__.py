__docformat__ = 'restructuredtext'

from mayaviViewer import MayaviViewer

__all__ = ["MayaviViewer"]

def make(*args, **kwargs):
    """
    A deprecated synonym for `MayaviViewer`
    """
    import warnings
    warnings.warn("'MayaviViewer' should be used instead of 'make'", DeprecationWarning, stacklevel=2)
    return MatplotlibViewer(*args, **kwargs)
