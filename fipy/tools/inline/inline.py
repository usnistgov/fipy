import sys

import Numeric 
import weave
from weave import converters

doingInline = 0

def doInline():
    global doingInline
    doingInline = 1
    
def dontInline():
    global doingInline
    doingInline = 0

def optionalInline(inlineFn, pythonFn, *args):
    global doingInline
    if doingInline:
	return inlineFn(*args)
    else:
	return pythonFn(*args)
	
def readInlineArgs():
    if "inline" in sys.argv:
	doInline()
    else:
	dontInline()

def runInline(code,**args):
    
    weave.inline(code,
		 args.keys(),
		 local_dict=args,
		 type_converters=converters.blitz,
		 compiler = 'gcc',
		 extra_compile_args =['-O3'])
			 
def runInlineLoop3(code_in,**args):
    
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
    
    return runInline(code, **args)
    
def runInlineLoop2(code_in,**args):
    
    code="""
    int i,j;
    for(i=0;i<ni;i++)
     {
      for(j=0;j<nj;j++)
       {
    """ + code_in + """
       }
     }
    """
    
    return runInline(code, **args)
    
def runInlineLoop1(code_in,**args):
    
    code="""
    int i;
    for(i=0;i<ni;i++)
     {
    """ + code_in + """
     }
    """
    
    print code
    
    return runInline(code, **args)

def runInlineLoop(code_in,**args):
    
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

