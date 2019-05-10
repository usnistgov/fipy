from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.terms.diffusionTerm import DiffusionTerm as ImplicitDiffusionTerm

__all__ = ["ImplicitDiffusionTerm"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]
