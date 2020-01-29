from __future__ import print_function
from __future__ import unicode_literals
from builtins import input
import os
from distutils.core import Command
from future.utils import text_to_native_str

from ._nativize import nativize_all

__all__ = [text_to_native_str("copy_script")]

class copy_script(Command):
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
    user_options = [nativize_all(u) for u in user_options]

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
                ans = input("The file '%s' already exists. Overwrite? [n] "%self.To)

            if ans is '':
                ans = 'no'

            if ("no".find(ans.lower()) is 0):
                self.To = input("Please give a name for the ouput file: ")
                self.finalize_options()

    def run(self):
        import imp
        import fipy.tests.doctestPlus

        mod = imp.load_source("copy_script_module", self.From)
        script = fipy.tests.doctestPlus._getScript(name = "copy_script_module")
        script = "\n\n## This script was derived from\n## '%s'\n\n%s"%(self.From, script)
        with open(self.To, "w") as f:
            f.write(script)

        print("Script code exported from '%s' to '%s'"%(self.From, self.To))
