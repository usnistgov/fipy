import sys

inlineCommentOn = False

if '--inline' in sys.argv[1:]:
    inlineFlagOn = True
    
    if '--verbose' in sys.argv[1:] or '-v' in sys.argv[1:]:
        inlineCommentOn = True
else:
    inlineFlagOn = False
