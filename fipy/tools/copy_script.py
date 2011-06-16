import os

class Copy_script():
    
    def __init__(self, To=None, From=None):
        self.To = To
        self.From = From

    def finalize_options(self):
        if self.From == None:
            raise "Please specify a '--From' input script file"
         
        if self.To == None:
            raise "Please specify a '--To' output script file"
            
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
        print 'script',script
        script = "#!/usr/bin/env python\n\n## This script was derived from\n## '%s'\n\n%s"%(self.From, script)
        raw_input('stopped 1')
        print 'script',script
        f = file(self.To, "w")
        f.write(script)
        f.close()
        
        print "Script code exported from '%s' to '%s'"%(self.From, self.To)
