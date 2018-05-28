#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "parser.py"
 #
 # Stolen from:
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import optparse
import sys, os

__all__ = ["parse"]

def parse(larg, action = None, type = None, default = None):
    """
    This is a wrapper function for the python `optparse` module.
    Unfortunately `optparse` does not allow command line arguments to
    be ignored. See the documentation for `optparse` for more
    details. Returns the argument value.

    :Parameters:
      - `larg`: The argument to be parsed.
      - `action`: `store` or `store_true` are possibilities
      - `type`: Type of the argument. `int` or `float` are possibilities.
      - `default`: Default value.

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
    elif 'FIPY_SOLVERS' in os.environ:
        return os.environ['FIPY_SOLVERS'].lower()
    else:
        return None
