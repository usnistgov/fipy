import sys

import Numeric 
##import weave

def _optionalInline(inlineFn, pythonFn, *args):
    if '--inline' in sys.argv[1:]:
	return inlineFn(*args)
    else:
	return pythonFn(*args)
	
def _runInline(code,**args):
    import weave

##    from weave.blitz_tools import blitz_type_factories	

##     print "in:", code
    
    weave.inline(code,
		 args.keys(),
		 local_dict=args,
##		 type_factories = blitz_type_factories,
##		 type_converters = weave_type_factories,
		 type_converters=weave.converters.blitz,
##		 type_converters=weave.weave_type_Factories
		 compiler = 'gcc',
		 verbose = 2,
		 extra_compile_args =['-O3'])
		 
##     print "out"
			 
def _runInlineLoop3(code_in,**args):
    
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
    
    return _runInline(code, **args)
    
def _runInlineLoop2(code_in,**args):
    
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
    
    return _runInline(code, **args)
    
def _runInlineLoop1(code_in,**args):
    
    code="""
    int i;
    for(i=0;i<ni;i++)
     {
    """ + code_in + """
     }
    """
    
    return _runInline(code, **args)

def _runInlineLoop(code_in,**args):
    
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
		 type_converters=weave.converters.blitz,
		 compiler = 'gcc',
		 extra_compile_args =['-O3'])

