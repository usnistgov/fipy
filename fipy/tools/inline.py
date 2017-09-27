__all__ = ["doInline"]

import inspect
import os
import sys

if '--inline' in [s.lower() for s in sys.argv[1:]]:
    doInline = True
else:
    doInline = 'FIPY_INLINE' in os.environ

_inlineFrameComment = 'FIPY_INLINE_COMMENT' in os.environ

def _getframeinfo(level, context=1):
    """
    Much faster alternative to `inspect.getouterframes(inspect.currentframe())[level]`
    """
    frame = inspect.currentframe()
    for l in range(level+1):
        frame = frame.f_back

    return (frame,) + inspect.getframeinfo(frame, context=context)

def _rawCodeComment(code, level=2):
    if _inlineFrameComment:
        finfo = _getframeinfo(level=level)

        # note:
        # don't use #line because it actually makes it harder
        # to find the offending code in both the C++ source and in the Python
        #line %d "%s"

        return '''
/*
    %s:%d

    %s
*/
        ''' % (finfo[1], finfo[2] - len(code.splitlines()), finfo[3])
    else:
        return ""

def _operatorVariableComment(canInline=True, level=3):
    if canInline and doInline and _inlineFrameComment:
        finfo = _getframeinfo(level=level)

        # note:
        # don't use #line because it actually makes it harder
        # to find the offending code in both the C++ source and in the Python
        #line %d "%s"

        if finfo[4] is not None:
            code = "\n".join(finfo[4])
        else:
            code = ""

        return '''
/*
%s:%d

%s
 */
        ''' % (finfo[1], finfo[2], code)
    else:
        return ""

def _runInline(code_in, converters=None, verbose=0, comment=None, **args):
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

    if comment is None:
        comment = _rawCodeComment(code_in)

    code = "\n" + comment + "\n" + code

    import weave

    for key in args.keys():
        if hasattr(args[key], 'dtype') and args[key].dtype.char == '?':
            args[key] = args[key].astype('B')

    weave.inline(code,
                 args.keys(),
                 local_dict=args,
                 type_converters=None, #weave.converters.blitz,
                 compiler = 'gcc',
                 force=0,
                 verbose = 0 or verbose,
                 extra_compile_args =['-O3'])

def _runIterateElementInline(code_in, converters=None, verbose=0, comment=None, **args):
    loops = """
int i;
for(i=0; i < ni; i++) {

"""

    enders = ""

    shape = args['shape']
    rank = len(shape) - 1
    for dim in range(rank):
        loops += "\t" * (dim + 1) + "for (vec[%(dim)d]=0; vec[%(dim)d] < shape[%(dim)d]; vec[%(dim)d]++) {\n" % {'dim': dim}
        enders += "\n" + "\t" * (rank - dim) + "}"

    enders += """

}
"""

    indent = "\t" * rank

    code = """
    #define ITEM(arr,i,vec) (arr[arrayIndex(arr##_array, i, vec)])

    int vec[%(rank)d];
    %(loops)s%(indent)s%(code_in)s%(enders)s

    #undef ITEM
    """ % locals()

    if comment is None:
        comment = _rawCodeComment(code_in)

    code = "\n" + comment + "\n" + code

    import weave

    for key in args.keys():
        if hasattr(args[key], 'dtype') and args[key].dtype.char == '?':
            args[key] = args[key].astype('B')

    weave.inline(code,
                 args.keys(),
                 local_dict=args,
                 type_converters=None, #weave.converters.blitz,
                 compiler = 'gcc',
                 force=0,
                 verbose = 0 or verbose,
                 extra_compile_args =['-O3'],
                 support_code="""

// returns the index (accounting for strides) of the tensor element vec
// in position i of array
//
// array holds a tensor at each position i
// vec identifies a particular element in that tensor
static int arrayIndex(PyArrayObject* array, int i, int vec[])
{
    int index = array->strides[array->nd-1] * i;

    if (vec != NULL) {
        int j;
        for (j=0; j < array->nd - 1; j++) {
            index += array->strides[j] * vec[j];
        }
    }

    return index / array->descr->elsize;
}
                 """)
