import sys
import Numeric 
import weave

def _optionalInline(inlineFn, pythonFn, *args):
    if '--inline' in sys.argv[1:]:
	return inlineFn(*args)
    else:
	return pythonFn(*args)
			 
def _runInline(code_in, converters=weave.converters.blitz, verbose=0, **args):
    #can add to dimList to increase dimensionality
    
    from fipy.variables.variable import Variable
    dimList = ['i', 'j', 'k']
    argsKeys = args.keys()
                
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

    
    #from scipy import weave

##    from weave.blitz_tools import blitz_type_factories	

##     print "in:", code
##    import weave
##
##    print 'code = ', code
##    print 'argDict =', args
##    for key in args.keys():
##        print 'key',key,'type(args[key])',type(args[key])
##        import MA

    weave.inline(code,
		 args.keys(),
      		 local_dict=args,
##		 type_factories = blitz_type_factories,
##		 type_converters = weave_type_factories,
##		 type_converters=None,
                 type_converters=converters,
##		 type_converters=weave.weave_type_Factories
		 compiler = 'gcc',
		 verbose = verbose,
		 extra_compile_args =['-O3'])
		 
