__all__ = ["PRINT"]

def PRINT(label, arg="", stall=True):
    import sys
    from fipy import parallel
    import time
    
    for procID in range(parallel.Nproc):
        if procID == parallel.procID:
            print >>sys.stderr, parallel.procID, label, arg
        sys.stderr.flush()
        if stall:
            time.sleep(0.1)
            parallel.Barrier()

    if stall:
        parallel.Barrier()
