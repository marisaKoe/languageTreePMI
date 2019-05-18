cimport numpy as np

cpdef np.float64_t nw(str x, str y, dict lodict, double gp1, double gp2)
cdef np.float64_t scoreNW(str x,str y, dict lodict, double gp1, double gp2)
cpdef np.ndarray[np.float64_t, ndim=2] compute_similarity(list l1List, list l2List, dict lodict, double gp1, double gp2)
cdef double gmean(np.ndarray[double, ndim=1] mylist)
cdef np.ndarray[np.float64_t, ndim=1] compute_ranks(np.ndarray[np.float64_t, ndim=1] diag, np.ndarray[np.float64_t, ndim=1] other)
cpdef tuple ldistNWPV(np.ndarray[np.float64_t, ndim=2] simMtr, double maxSim, double minSim)