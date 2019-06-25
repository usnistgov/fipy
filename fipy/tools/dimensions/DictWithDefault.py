from __future__ import unicode_literals
# A dictionary with default values for non-existing entries

__all__ = []

import copy

class _DictWithDefault(dict):
    """Dictionary with default values

    Instances of this class act like standard Python dictionaries,
    except that they return a *copy* of |default| for a key that
    has no associated value.
    """

    def __init__(self, default):
        self.default = default

    def __getitem__(self, key):
        try:
            item = dict.__getitem__(self, key)
        except KeyError:
            item = copy.copy(self.default)
            self[key] = item
        return item

    def __delitem__(self, key):
        try:
            del self.data[key]
        except KeyError:
            pass
