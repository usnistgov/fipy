import Numeric
import weave
from weave import converters

doingInline = 0

def readInlineArgs(argList):
    if "inline" in argList:
	doInline()
    else:
	dontInline()

def doInline():
    doingInline = 1
    
def dontInline():
    doingInline = 0

def optionalInline(inlineFn, pythonFn, *args):
    if doingInline:
	return inlineFn(*args)
    else:
	return pythonFn(*args)
    
def runInline(code_in,**args):
    
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
    
    weave.inline(code,
                 args.keys(),
                 local_dict=args,
                 type_converters=converters.blitz,
                 compiler = 'gcc',
                 extra_compile_args =['-O3'])

