import numpy as np
cimport numpy as np

cdef extern from "math.h":
    double log(double)

cdef extern from "math.h":
    double exp(double)

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


cpdef np.float64_t nw(str x, str y, dict lodict, double gp1, double gp2):
    """
    Needleman-Wunsch algorithm for pairwise string alignment
    with affine gap penalties.
    'lodict' must be a dictionary with all symbol pairs as keys
    and match scores as values.
    gp1 and gp2 are gap penalties for opening/extending a gap.
    Returns the alignment score and one optimal alignment.
    """
    #length of the words
    cdef int n,m
    cdef list pointers
    cdef np.ndarray[np.float64_t, ndim=2] dp
    cdef int i,j
    cdef np.float64_t match, insert, delet, temp
    n = len(x)
    m = len(y)

    #dp = [[0.0 for i in range(m+1)] for j in range(n+1)]
    dp = np.zeros((n+1,m+1), DTYPE)
    pointers = [[0 for i in range(m + 1)] for j in range(n + 1)]

    dp[1][0] = gp1
    pointers[1][0] = 1
    dp[0][1] = gp1
    pointers[0][1] = 2
    for i in range(2,n+1):
        dp[i][0] = dp[i-1][0]+gp2
        pointers[i][0] = 1
    for j in range(2,m+1):
        dp[0][j] = dp[0][j-1]+gp2
        pointers[0][j] = 2
    for i in range(1, n+1, 1):
        for j in range(1, m+1, 1):

            match = dp[i-1][j-1]+lodict[(x[i-1],y[j-1])]

            temp = gp1
            if pointers[i-1][j] == 1:
                temp = gp2
            insert = dp[i-1][j]+temp

            temp = gp1
            if pointers[i][j-1] == 2:
                temp = gp2
            delet = dp[i][j-1]+temp

            if match >= insert and match >= delet:
                dp[i][j] = match
                pointers[i][j] = 0
            elif insert >= match and insert >= delet:
                dp[i][j] = insert
                pointers[i][j] = 1
            else:
                dp[i][j] = delet
                pointers[i][j] = 2

    return dp[-1][-1]

cdef np.float64_t scoreNW(str x,str y, dict lodict, double gp1, double gp2):
    """
    Calculate NW score for two words, possibly synonyms
    :param x: word 1
    :param y: word 2 
    :param lodict: dictionary of similarity scores
    :param gp1: gap opening penalty
    :param gp2: gap extension penalty
    :return: similarity score
    """
    cdef str  w
    cdef list x1,y1

    if x == '0' or y == '0':
        # missing entry
        return np.nan

    # synonyms are marked with - in the asjp matrix, split the synonyms,
    # align all combinations and get the maximal similarity score
    x1=[w for w in x.split('-') if not '%' in w]
    y1=[w for w in y.split('-') if not '%' in w]

    return max([nw(xx,yy,lodict,gp1,gp2) for xx in x1 for yy in y1])

cpdef np.ndarray[np.float64_t, ndim=2] compute_similarity(list l1List, list l2List, dict lodict, double gp1, double gp2):
    """
    Computes a similarity matrix of two word lists using Needleman-Wunsch
    :param l1List: word list of language 1
    :param l2List: word lsit of language 2
    :param lodict: dictionary of similarity scores
    :param gp1: gap opening penalty
    :param gp2: gap extension penalty 
    :return: similiarity matrix
    """
    cdef Py_ssize_t x, y
    cdef np.ndarray[np.float64_t, ndim=2] simMtr = np.empty([len(l2List), len(l1List)], DTYPE)


    for y in range(len(l2List)):
        for x in range(len(l1List)):
            simMtr[y,x] = scoreNW(l1List[x],l2List[y], lodict, gp1, gp2)

    return simMtr

cdef double gmean(np.ndarray[double, ndim=1] mylist):
    """
    calculate geometric mean
    :param mylist: list to calculate geometric mean
    :return: geometric mean
    """
    cdef double product = 0.0
    cdef Py_ssize_t i

    if len(mylist) == 1:
        return mylist[0]

    for i in range(len(mylist)):
        product += log(mylist[i])

    return exp((1./len(mylist))*product)


cdef np.ndarray[np.float64_t, ndim=1] compute_ranks(np.ndarray[np.float64_t, ndim=1] diag, np.ndarray[np.float64_t, ndim=1] other):
    """
    compute ranks of diagonal scores
    :param diag: diagonal alignment scores
    :param other: off-diagonal alignment scores
    :return: ranks of diagonal scores
    """
    cdef np.ndarray[np.float64_t, ndim=1] ranks = np.empty(len(diag), DTYPE)
    cdef Py_ssize_t i
    cdef double x
    cdef int geq, g, n
    cdef np.float64_t value


    other = np.sort(other, kind="quicksort")
    n = len(other)
    for i in range(len(diag)):
        x = diag[i]

        # greater and equal
        geq = n - np.searchsorted(other, x)

        # greater
        g = n - np.searchsorted(other, x, "right")

        value = gmean(1.0+np.arange(g, 1.0+geq))
        ranks[i] = value

    return ranks

cpdef tuple ldistNWPV(np.ndarray[np.float64_t, ndim=2] simMtr,
                       double maxSim, double minSim):
    """
    distNWPV is the distance measure that is called dERC/PMI in Jaeger (2013)
    :param simMtr: similarity matrix (calculated from compute_similarity)
    :param maxSim: maximal similarity score
    :param minSim: minimal similarity score
    :return: distance measure between languages
    """
    cdef np.ndarray[np.float64_t, ndim=1] dg, cmpr, ranks
    cdef double stc, sim


    dg = np.diag(simMtr) #get diagonal
    dg = dg[np.isnan(dg)==False] #get non nan subset
    np.fill_diagonal(simMtr,np.nan) #set diagonal to nan
    cmpr = simMtr[np.isnan(simMtr)==False] # flatten to vector
    ranks = compute_ranks(dg, cmpr) # get ranks
    stc = np.mean(-np.log(ranks/(1+len(cmpr)))) # mean of log ranks
    sim = (stc-1)*np.sqrt(len(dg)) # similarity score

    return (maxSim-sim)/(maxSim-minSim), ranks
