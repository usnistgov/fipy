import Numeric
import weave
from weave import converters

class Inline:

    doingInline = 0
    
    def doInline(self):
        print 'in doInline'
        self.doingInline = 1
    
    def dontInline(self):
        print 'in dontInline'
        self.doingInline = 0

    def readInlineArgs(self,argList):
        if "inline" in argList:
            self.doInline()
        else:
            self.dontInline()

    def optionalInline(self,inlineFn, pythonFn, *args):
        if self.doingInline:
            return inlineFn(*args)
        else:
            return pythonFn(*args)
    

    def runInline(code_in, **args):

        code="""
        int i,j,k;
        for(i=0;i<ni;i++)
         {
          for(j=0;j<nj;j++)
           {
            for(k=0;k<nk;k++)
             {
        """ + code_in + """
             }
            }
           }
        """

        weave.inline(
            code,
            args.keys(),
            local_dict = args,
            type_converters = converters.blitz,
            compiler = 'gcc',
        extra_compile_args = ['-O3'])

inl = Inline()
