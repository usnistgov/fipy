def PRINT(label, *args, **kwargs):
    stall = kwargs.get('stall', True)
    
    import sys
    from fipy import parallel
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    import time
    
    for procID in range(parallel.Nproc):
        if procID == parallel.procID:
            print >>sys.stderr, parallel.procID, label, args
        sys.stderr.flush()
        if stall:
            time.sleep(0.1)
            comm.Barrier()

    if stall:
        comm.Barrier()
