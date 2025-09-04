from builtins import range
__all__ = ["PRINT"]

def PRINT(label, arg="", stall=True):
    """Display `label` and `arg` on each MPI rank.

    Annotate with rank number.
    If `stall` is true, ensure all ranks print before proceeding.
    """
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
