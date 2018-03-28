## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "copy_script.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: Andrew Acquaviva <andrewa@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
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

import os
from distutils.core import Command

__all__ = ["Copy_script"]

class Copy_script(Command):
    description = "copy an example script into a new editable file"

    # List of option tuples: long name, short name (None if no short
    # name), and help string.
    user_options = [
        # Select installation scheme and set base director(y|ies)
        ('From=', None,
         "path and file name containing script to copy"),
        ('To=', None,
         "path and file name to save script to")
     ]

    def initialize_options(self):
        self.From = None
        self.To = None

    def finalize_options(self):
        if self.From == None:
            raise SyntaxError("Please specify a '--From' input script file")

        if self.To == None:
            raise SyntaxError("Please specify a '--To' output script file")

        if os.path.exists(os.path.expanduser(self.To)):
            ans = "junk"

            while (len(ans) > 0) and ("yes".find(ans.lower()) is not 0) and ("no".find(ans.lower()) is not 0):
                ans = raw_input("The file '%s' already exists. Overwrite? [n] "%self.To)

            if ans is '':
                ans = 'no'

            if ("no".find(ans.lower()) is 0):
                self.To = raw_input("Please give a name for the ouput file: ")
                self.finalize_options()

    def run(self):
        import imp
        import fipy.tests.doctestPlus

        mod = imp.load_source("copy_script_module", self.From)
        script = fipy.tests.doctestPlus._getScript(name = "copy_script_module")
        script = "#!/usr/bin/env python\n\n## This script was derived from\n## '%s'\n\n%s"%(self.From, script)
        f = file(self.To, "w")
        f.write(script)
        f.close()

        print "Script code exported from '%s' to '%s'"%(self.From, self.To)
