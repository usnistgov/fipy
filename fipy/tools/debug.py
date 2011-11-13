__all__ = ["PRINT"]

def PRINT(label, arg="", stall=True):
    import sys
    from fipy import parallel
    import time
    
    from six import print_
    
    for procID in range(parallel.Nproc):
        if procID == parallel.procID:
            print_(parallel.procID, label, arg, file=sys.stderr)
        sys.stderr.flush()
        if stall:
            time.sleep(0.1)
            parallel.Barrier()

    if stall:
        parallel.Barrier()
