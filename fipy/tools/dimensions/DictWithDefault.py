# A dictionary with default values for non-existing entries

import UserDict, copy

class DictWithDefault(UserDict.UserDict):

    """Dictionary with default values

    Constructor:  DictWithDefault(|default|)

    Instances of this class act like standard Python dictionaries,
    except that they return a *copy* of |default| for a key that
    has no associated value.
    """

    def __init__(self, default):
	self.data = {}
	self.default = default

    def __getitem__(self, key):
	try:
	    item = self.data[key]
	except KeyError:
	    item = copy.copy(self.default)
	    self.data[key] = item
	return item

    def __delitem__(self, key):
	try:
	    del self.data[key]
	except KeyError:
	    pass
