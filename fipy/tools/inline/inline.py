import Numeric
import weave
from weave import converters

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

