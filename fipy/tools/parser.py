from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

import optparse
import sys, os

__all__ = ["parse"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

def parse(larg, action = None, type = None, default = None):
    """
    This is a wrapper function for the python `optparse` module.
    Unfortunately `optparse` does not allow command line arguments to
    be ignored. See the documentation for `optparse` for more
    details. Returns the argument value.

    Parameters
    ----------
    larg : str
        Argument to be parsed.
    action : {`'store'`, `'store_true'`, `'store_false'`, `'store_const'`, `'append'`, `'count'`, `'callback'`}
        Basic type of action to be taken when this argument is encountered at the command line.
        See https://docs.python.org/2/library/argparse.html#action
    type : type
        Type to which the command-line argument should be converted
    default
        Value produced if the argument is absent from the command line
    """
    sarg = None
    tmpparser = optparse.OptionParser(option_list = [
        optparse.make_option(sarg, larg, action = action, type = type, dest = 'dest', default = default)],
                                      conflict_handler = 'resolve')

##    optparse.make_option('-e', '--numberOfElements', action = 'store', type = 'int', dest = 'Nele', default = numberOfElements),
##    optparse.make_option('-n', '--numberOfSteps', action = 'store', type = 'int', dest = 'steps', default = numberOfSteps),
##    optparse.make_option('-p', '--quiet', action = 'store_true', dest = 'quiet', default = False),
##    optparse.make_option('-i', '--inline', action = 'store_true', dest = 'inline', default = False)])

    sysargs = []
    for arg in sys.argv:
        if larg in arg:
            sysargs.append(arg)

    (options, args) = tmpparser.parse_args(sysargs)

    return options.dest

def _parseSolver():
    args = [s.lower() for s in sys.argv[1:]]
    # any command-line specified solver takes precedence over environment variables
    if '--no-pysparse' in args:
        return "no-pysparse"
    elif '--trilinos' in args:
        return "trilinos"
    elif '--pysparse' in args:
        return "pysparse"
    elif '--pyamg' in args:
        return 'pyamg'
    elif '--scipy' in args:
        return 'scipy'
    elif '--pyamgx' in args:
        return 'pyamgx'
    elif '--petsc' in args:
        return 'petsc'
    elif 'FIPY_SOLVERS' in os.environ:
        return os.environ['FIPY_SOLVERS'].lower()
    else:
        return None
