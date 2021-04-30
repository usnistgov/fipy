import logging
import numpy as np
from scipy.sparse.linalg import splu
from scipy.io import mmread

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s:%(levelname)s:%(message)s')

def L2norm(vec):
    return np.sqrt(np.sum(vec**2))

def error(L, x, b):
    return L * x - b

def solve(var):
    logging.debug("{}:START".format(var))

    L = mmread("{}-matrix.mtx".format(var))
    x = np.loadtxt("{}-xvec.dat".format(var))
    b = np.loadtxt("{}-bvec.dat".format(var))

    logging.debug("{}:LOAD".format(var))

    diag = L.diagonal()
    maxdiag = max(np.absolute(diag))

    L = L * (1 / maxdiag)
    b = b * (1 / maxdiag)

    logging.debug("{}:NORMALIZE".format(var))

    LU = splu(L.asformat("csc"),
              diag_pivot_thresh=1.,
              relax=1,
              panel_size=10,
              permc_spec="COLAMD",
              options=dict(PrintStat=True))

    logging.debug("{}:SETUP".format(var))

    error0 = L2norm(error(L, x, b))

    for iteration in range(10):
        errorVector = error(L, x, b)
        if (L2norm(errorVector) / error0) <= 1e-10:
            break

        xError = LU.solve(errorVector)
        x[:] = x - xError

        logging.debug("{}:SOLVE:{}".format(var, iteration))

    print("var: {}".format(var))
    print("iterations:", iteration+1)
    print("residual:", L2norm(errorVector))
    print("x", "min", min(abs(x)), "max", max(abs(x)))

    logging.debug("{}:END".format(var))

# solve("theta")
solve("phase")
# solve("heat")
