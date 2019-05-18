'''
Created on 15.12.2017

@author: marisa

TODO: set variable in sample_wrdLst


'''
import codecs, timeit
import numpy as np

from itertools import islice
from collections import defaultdict
from multiprocessing import Pool

import cython_code.code.nw as c

##angepasst auf NELex, values defined in Jaeger 2013
#minSim = -np.sqrt(1016)
#maxSim = (np.log(1016*(1015)+1)-1)*np.sqrt(1016)

##200 concpets
minSim = -np.sqrt(200)
maxSim = (np.log(200*(199)+1)-1)*np.sqrt(200)

def create_bootstrapMtx(simMtr):
    '''
    Creates 1000 samples of the wordlists in nelex.
    The wordlists are sampled with replacement and according to the concepts (same concepts in the lists of samples).
    Those samples can be used for the traditional bootstrap method.
    :return bootstrap_dict: dict with key=number of sample value=dict with key=language (iso-code) value=sample word lists with replacement
    '''
    
    ##get all the language pairs in a list
    lang_pair_list = simMtr.keys()
    number_concepts = len(simMtr.values()[0])
    #print number_concepts
    ##get the number of concepts to create the indice list
    #number_concpets = len(nelex_dict.values()[0])
    
    ##initialize the overall bootstrap dict
    bootstrap_dict = defaultdict()
    n = 1000
    for number_sample in range(n):
        start = timeit.default_timer()
        ##creates a list with indices between the range of 1016, the list is 1016 chars long -> do this 1000 times
        idx_list = np.random.choice(np.arange(number_concepts), number_concepts, replace=True)
        #idx_list = np.random.choice(np.arange(10), 10, replace=True)
        ##initialize the distance dict for each sample
        sample_dist_dict = defaultdict(lambda: defaultdict(float))
        #print "index list run", number_sample, idx_list
        ##for each pair in the list of language pairs
        for pair in lang_pair_list:
            ##create the matrix
            matrix=simMtr[pair]
            
            ##initialize the bootstrap array in numpy 
            bootstrap_array = np.zeros((len(idx_list),len(idx_list)))
             
            ##for the index (in the index list) and index of the word in the simMatrix
            for idx_idxList, idx_word in enumerate(idx_list):
                ##for each number in the index list, get the index and the index (number) of the word
                for jdx_idxList, jdx_word in enumerate(idx_list):
                    
                    ##create the bootstrap array out of the simMaitrx with all indices and set it in the matrix
                    ##we dont have to compute the similarities between the words again, and therefore the program runs faster
                    bootstrap_array[idx_idxList][jdx_idxList] = matrix[idx_word][jdx_word]
            

            
            #print bootstrap_array
            ##get the distance between a pair of languages and set it in the distance matrix
            dist_score, ranks = c.ldistNWPV(bootstrap_array, maxSim, minSim)
            #print "dist score", dist_score
            #stop = timeit.default_timer()
            #print stop-start
            sample_dist_dict[pair[0]][pair[1]] = dist_score
            sample_dist_dict[pair[1]][pair[0]] = dist_score
             
         
        ##append the sample dictionary to the bootstrap dictionary
        bootstrap_dict[number_sample] = sample_dist_dict
        stop = timeit.default_timer()
        print stop-start, " still computing bootstrap matrices"
     
    ##write all distance matrices in a folder
    out_path_mtx = "output/bootstrap/phylip/nelexDstMtx"
    ##write all matrices at once in seperate files into a folder
    for n_sample, dstMtx in bootstrap_dict.items():
        fout1 = codecs.open(out_path_mtx+"+"+str(n_sample)+".phy","wb","utf-8")
        writePhy(dstMtx, fout1)
        


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
    f = open("input/simMtx.phy")
    raw_data = f.readlines()
    f.close()
    
    simMtx = defaultdict()
    for l in raw_data:
        line = l.split("\t")
        langPair = line[0]
        simArray = line[1:]
        lp = tuple(langPair.split(","))
        
        newSimArray = list()
        for s in simArray:
            s1 = s.split()
            #print "s1", s1
            newList = [float(i) for i in s1]
            #print "newList", newList
            if len(newList) !=0:
                newSimArray.append(newList)
        simAr = np.array(newSimArray)
        simMtx[lp] = simAr
    
    
    
    create_bootstrapMtx(simMtx)
    
    
    
    
    
    #create_bootstrapMtx()
    
    
    