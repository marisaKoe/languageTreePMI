'''
Created on 15.12.2017

@author: marisa

'''

import os, codecs, timeit

import itertools
import numpy as np
from Bio import pairwise2
from scipy import stats, special
from collections import defaultdict
from itertools import combinations

import cython_code.code.nw as c


##################global variables for the algorithm#####################

##gap penalties for the Needleman-Wunsch algorithm
gp1 = -2.49302792222
gp2 = -1.70573165621
 

##angepasst auf NELex, values defined in Jaeger 2013
#minSim = -np.sqrt(1016)
#maxSim = (np.log(1016*(1015)+1)-1)*np.sqrt(1016)

##200 concpets
minSim = -np.sqrt(200)
maxSim = (np.log(200*(199)+1)-1)*np.sqrt(200)

def compute_dm(nelex_dict, lodict):
    '''
    computes the distance matrix and writes it in a file
    
    :param nelex_dict:the dictionary for the database key=language (iso-code) value=list of words
    :param out_path:the path for the distance matrix
    '''
    
    
    ##get a list with all languages
    list_langs = nelex_dict.keys()
    
    ##initialize the dictionary for the distances
    dist_dict = defaultdict(lambda: defaultdict(float))
    ##initialize overall simMtr for all language pairs and their word lists (important for sampling afterwards)
    overall_simMtr = defaultdict()
    ##create language pairs
    count = 0
    for lang_pair in combinations(list_langs,r=2):
        start = timeit.default_timer()
        ##get the languages
        l1,l2 = lang_pair
        count += 1
        ##get the word lists for each language
        l1List = nelex_dict[l1]
        l2List = nelex_dict[l2]

        ##compute the similarity matrix between the wordlists
        simMtr = c.compute_similarity(l1List, l2List, lodict, gp1, gp2)
        #simMtr = compute_similarity_matrix(l1List, l2List, lodict, gp1, gp2)
        
        overall_simMtr[lang_pair]=simMtr
        simMtr_2 = np.copy(simMtr)
        ##computes the distance between the languages
        dist_score, ranks = c.ldistNWPV(simMtr_2, maxSim, minSim)
      
        ##fill the dictionary with the values
        dist_dict[l1][l2] = dist_score 
        dist_dict[l2][l1]=dist_score
        
        stop = timeit.default_timer()
        print stop-start, " still computing matrices", count
    
    
    return dist_dict, overall_simMtr
    



# #############################distance methods###################################

def sscore(a,b, lodict, go1=gp1, gp2=gp2):
    return pairwise2.align.globalds(a,b,lodict,gp1,gp2)[0][2]


def nw(x,y,lodict,gp1=gp1,gp2=gp2):
    """
    Needleman-Wunsch algorithm for pairwise string alignment
    with affine gap penalties.
    'lodict' must be a dictionary with all symbol pairs as keys
    and match scores as values.
    gp1 and gp2 are gap penalties for opening/extending a gap.
    Returns the alignment score and one optimal alignment.
    """
    #length of the words
    n,m = len(x),len(y)
  
    dp = [[0.0 for i in range(m+1)] for j in range(n+1)]
    pointers = [[0 for i in range(m + 1)] for j in range(n + 1)]
  
    dp[1][0] = gp1
    pointers[1][0] = 1
    dp[0][1] = gp1
    pointers[0][1] = 2
    for i in xrange(2,n+1):
        dp[i][0] = dp[i-1][0]+gp2
        pointers[i][0] = 1
    for j in xrange(2,m+1):
        dp[0][j] = dp[0][j-1]+gp2
        pointers[0][j] = 2
    for i, j in itertools.product(range(n+1)[1:], range(m+1)[1:]):
  
        match = (dp[i-1][j-1]+lodict[x[i-1],y[j-1]], 0)
        insert = (dp[i-1][j]+(gp2 if pointers[i-1][j]==1 else gp1), 1)
        delet = (dp[i][j-1]+(gp2 if pointers[i][j-1]==2 else gp1), 2)
        dp[i][j],pointers[i][j] = max([match,insert,delet], key=lambda item:item[0])
  
    return dp[-1][-1]
   
   
   
def scoreNW(x,y,lodict,gp1=gp1,gp2=gp2):
    '''
    Function to handle synonmys in the data.
    Compute the Needleman-Wunsch score (similarity) between all combination of words in the lists and take the pair with the maximal similarity.
    This is important for synonyms, if a language only has one word for the concept, this score is taken.
    :param x:word(s) for the concept for one language
    :param y:word(s) for the concept for the other language
    :param lodict:the dictionary with the sounds and pmis
    :param gp1:gap penalty (opening gap)
    :param gp2:gap penalty (extending gap)
    '''
    if '0' in [x,y]: return np.nan
    ##synonyms are marked with - in the asjp matrix, split the synonyms, align all combinations and get the maximal similarity score
    x1=[w for w in x.split('-') if not '%' in w]
    y1=[w for w in y.split('-') if not '%' in w]
    return max([nw(xx,yy,lodict,gp1,gp2) for xx in x1 for yy in y1])
   
   
# distNWPV is the distance measure that is called dERC/PMI in Jaeger (2013)
def compute_similarity_matrix(l1List,l2List,lodict,gp1=gp1,gp2=gp2):
    '''
    computes the similarity matrix between all words in the list of languages
    :param l1List:word list for language 1
    :param l2List:word list for language 2
    :param lodict:dictionary with logodds for the sounds
    :param gp1:gap penalty (opening gap)
    :param gp2:gap penalty (extending gap)
    :return simMtr: similarity matrix for the language lists
    '''
    ##give each word in the lists to the NW function and create a similarity matrix
    sim_matrix = np.array([[scoreNW(x,y,lodict,gp1=gp1,gp2=gp2) for x in l1List] for y in l2List])
      
    return sim_matrix
      
  
# distNWPV is the distance measure that is called dERC/PMI in Jaeger (2013)
def compute_ldistNWPV(simMtr):
    '''
    together with compute_similarity_matrix => ldistNWPV is the distance measure that is called dERC/PMI in Jaeger (2013)
    dERC = distance based on corrected evidence of relatedness
    based on PMI scores, which is better for the computation (evidence in Jaeger 2013)
    :param simMtr:similarity matrix from compute_similarity_matrix
    :return: distance between the languages
    '''
    ##get the diagonal of the matrix
    dg = np.diag(simMtr)
    ##tests element-wise if there are non nan, get the subset
    dg = dg[np.isnan(dg)==False]
    ##sets diagonal values to nan
    np.fill_diagonal(simMtr,np.nan)
    #print "filled diagonals"
    ##compromise the matrix to a list/vector
    cmpr = simMtr[np.isnan(simMtr)==False]
      
    ##gives for each word in the diagonal (synonyms) a rank
    ranks = np.array([stats.gmean(1.+np.arange(sum(cmpr>x),1.+sum(cmpr>=x))) for x in dg],np.double)
    stc = np.mean(-np.log(ranks/(1+len(cmpr))))
    sim = (stc-1)*np.sqrt(len(dg))
    return (maxSim-sim)/(maxSim-minSim)

###############writing################################


def writePhy(dm,f):
    '''
    Writes the distance matrix into a file of format phylip
    :param dm: the dictionary with the distances
    :param f: the output file name
    '''
    ntaxa = len(dm)
    f.write(str(ntaxa)+"\n")
        
    i=0
    for ka in dm:
        row=ka+" "
        for kb in dm:
            row=row+" "+str(dm[ka][kb])
        f.write(row+"\n")
    f.close()






if __name__ == '__main__':
    
    
    compute_dm()
    

    
    
    
    
    
    
    
    