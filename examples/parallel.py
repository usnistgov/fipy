from fipy import parallel, Grid1D
mesh = Grid1D(nx=10)
print "%d cells on processor %d of %d" \
  % (mesh.getNumberOfCells(), parallel.procID, parallel.Nproc)
