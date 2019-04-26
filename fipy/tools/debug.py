from __future__ import print_function
__all__ = ["PRINT"]

def PRINT(label, arg="", stall=True):
    import sys
    from fipy import parallelComm
    import time

    for procID in range(parallelComm.Nproc):
        if procID == parallelComm.procID:
            print(parallelComm.procID, label, arg, file=sys.stderr)
        sys.stderr.flush()
        if stall:
            time.sleep(0.1)
            parallelComm.Barrier()

    if stall:
        parallelComm.Barrier()
