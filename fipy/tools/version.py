"""Shim for version checking

`distutils.version` is deprecated, but `packaging.version` is unavailable
in Python 2.7
"""
from __future__ import unicode_literals

__docformat__ = 'restructuredtext'


__all__ = ["Version", "parse_version"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

try:
    from packaging.version import Version, parse as parse_version
except ImportError:
    from distutils.version import StrictVersion as Version, StrictVersion as parse_version    
