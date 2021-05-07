import logging
from scipy.sparse.linalg import splu
from scipy.io import mmread

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s:%(levelname)s:%(message)s')

logging.debug("START")

L = mmread("phase-matrix.mtx")

logging.debug("mmread")

Lcsc = L.asformat("csc")

logging.debug("asformat_csc")

LU = splu(Lcsc,
          diag_pivot_thresh=1.,
          relax=1,
          panel_size=10,
          permc_spec="COLAMD")

logging.debug("splu")

logging.debug("END")
