import sys

if '--inline' in sys.argv[1:]:
    inlineFlagOn = True
else:
    import os

    if os.environ.has_key('FIPY_INLINE'):
        if os.environ['FIPY_INLINE'].lower() == 'true':
            inlineFlagOn = True
        elif os.environ['FIPY_INLINE'].lower() == 'false':
            inlineFlagOn = False
        else:
            raise ImportError, 'Unknown setting for FIPY_INLINE: %s, should be true or false' % os.environ['FIPY_INLINE']
    else:
        inlineFlagOn = False
