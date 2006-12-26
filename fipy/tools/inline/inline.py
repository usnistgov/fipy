import sys
from fipy.tools import numerix 

def _optionalInline(inlineFn, pythonFn, *args):
    if '--inline' in sys.argv[1:]:
	return inlineFn(*args)
    else:
	return pythonFn(*args)
			 
def _runInline(code_in, converters=None, verbose=0, **args):
    argsKeys = args.keys()
    dimList = ['i', 'j', 'k']
          
    if 'ni' in argsKeys:
        dimensions = 1
        if 'nj' in argsKeys:
            dimensions = 2
            if 'nk' in argsKeys:
                dimensions = 3
    else:
        dimensions = 0
    
    if dimensions == 0:
        code = """ { %s } """ % code_in 
    else:
        loops = """"""
        enders = """"""
        declarations = []
        for dim in range(dimensions):
            d = dimList[dim]
            declarations.append(d)
            loops += "\t" * dim + "for(%s=0;%s<n%s;%s++) {\n" % (d,d,d,d)
            enders += "\n" + "\t" * (dimensions - dim -1) + "}"
        code = 'int ' + ','.join(declarations) + ';\n' + loops + "\t" * dimensions + code_in + enders

    from scipy import weave
        
    weave.inline(code,
                 args.keys(),
                 local_dict=args,
                 type_converters=converters or weave.converters.blitz,
                 compiler = 'gcc',
                 verbose = verbose,
                 extra_compile_args =['-O3'])
