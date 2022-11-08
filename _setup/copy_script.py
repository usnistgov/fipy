from __future__ import print_function
import os
from distutils.core import Command

__all__ = ["copy_script"]

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
                from builtins import input
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
