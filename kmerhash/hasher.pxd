#cython: language_level=3

import numpy as np
cimport numpy as np


cpdef np.ndarray[np.uint64_t, ndim=1] kmerhasher(str seq, int kmer_length)


