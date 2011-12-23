"""The bibstyles package provides classes and styles for formatted text output.

This is the __init__.py for bibstyles package
"""

from . import default
from . import example_numbered
from . import jasss_style

_named = dict(
    default = default,
    example_numbered = example_numbered,
    jasss_style = jasss_style,
    jasss = jasss_style)

def from_name(name):
    return _named[name]
