import numpy as np
from scipy.sparse.linalg import splu
from scipy.io import mmread

def L2norm(vec):
    return np.sqrt(np.sum(vec**2))

def error(L, x, b):
    return L * x - b

L = mmread("matrix.mtx")
x = np.loadtxt("xvec.dat")
b = np.loadtxt("bvec.dat")

diag = L.diagonal()
maxdiag = max(np.absolute(diag))

L = L * (1 / maxdiag)
b = b * (1 / maxdiag)

LU = splu(L.asformat("csc"),
          diag_pivot_thresh=1.,
          relax=1,
          panel_size=10,
          permc_spec="COLAMD",
          options=dict(PrintStat=True))

error0 = L2norm(error(L, x, b))

for iteration in range(10):
    errorVector = error(L, x, b)
    if (L2norm(errorVector) / error0) <= 1e-10:
        break

    xError = LU.solve(errorVector)
    x[:] = x - xError

print("iterations:", iteration+1)
print("residual:", L2norm(errorVector))
print("x", "min", min(abs(x)), "max", max(abs(x)))