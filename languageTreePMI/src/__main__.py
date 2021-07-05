'''
Created on 15.12.2017

@author: marisa

Method read_original_nelex:
------------------------------
Write method to read in original NELex format and create an ASJP matrix format out of it. 
+ synonyms are written together and separated with -
+ matrix needs to be computed with the new version, there are mistakes in the old one
'''
import codecs, timeit 
import itertools as it
import numpy as np


from collections import defaultdict
import compute_distmtx, reconstruct_trees, traditional_bootstrap, create_consensus_paup




def read_nelex():
    '''
    read the nelex database in the asjp matrix format from Gerhard (downloaded from trac/gjaeger
    :return nelex_dict: dictionary with key=languages(iso-code) value = list of words
    '''
    
    ##read the matrix form Gerhard for the nelex sample in asjp format
    f = open("input/nelexAsjpMtx.txt")
    raw_data = f.readlines()
    f.close()
    
    
    ##dictionary with key = language/iso-code value=list with all asjp_words
    nelex_dict = defaultdict()
    ##for each line in the matrix, get the language and the list of words and append it to a dictionary
    for l in raw_data:
        l = l.strip().split("\t")
        lang = l[0]
        word_list = l[1:]
        
        nelex_dict[lang]=word_list
    
    return nelex_dict


def read_nelex_original():
    '''
    read the nelex database, which can be downloaded
    prepares the data:
    + include missing entries with a placeholder 0 (the strings should be of same length for the alignments)
    + join synonyms with an eliminater: -
    :return nelex_dict: dictionary with key= languages (iso-code) and value=list of words for this language
    '''
    
    f = open("input/northeuralex-cldf.tsv")
    raw_data = f.readlines()
    f.close()
    
    ##read only the IE languages to get a smaller sample
    with open("input/nelex-ie-langs.txt") as f1:
        ie_langs = f1.read().splitlines()
        
    ##read the top 200 concepts of nelex to get a smaller sample for the bootstrap
    top_concepts = []
    with open("input/Top200Nelex.csv") as f2:
        raw_concepts = f2.readlines()
        for line in raw_concepts[1:]:
            concept = line.strip().split("\t")[0]
            if " " in concept:
                concept = concept.replace(" ","")
            top_concepts.append(concept)
        
    
    
    ##dictionary with key = languages value=dict with key=concept value= list of words
    nelex_original_dict = defaultdict(lambda: defaultdict(list))
    unique_concepts = list()
    for l in raw_data[1:]:
        line = l.split("\t")
        ##glottocode
        #glot = line[0]
        ##iso code
        iso_code = line[1]
        
        ##concept
        concept = line[2]
        if " " in concept:
            concept = concept.replace(" ","")
        ##ipa
        #ipa = line[3]
        ##asjp
        asjp_word = line[4]
        ##sca
        #sca = line[5]
        ##dolgo
        #dolgo = line[6]
        
        ##make a list with all unique concepts
        if not concept in unique_concepts:
            unique_concepts.append(concept)
        
        ##create a dicitonary with the original data, NOTE: missing entries are just left out, need to be inserted afterwards for the specific concept
        if iso_code in ie_langs:
#             ##!!!!!!!!!!!!!!!!!!!comment out if the full concepts are used
#             if concept in top_concepts:
            nelex_original_dict[iso_code][concept].append(asjp_word)
     
    
    ##new dict in form of asjp matrix format form Gerhard   
    nelex_dict = defaultdict()
    
    ##for each language, concept dictionary
    for lang, con_dict in nelex_original_dict.items():
        #print len(con_dict)
        word_list_lang = list()
        '''comment in and out depending on the number of concepts used for the experiment'''
        ##go through the list of unique concepts if the 1016  concepts are used
        #for concept in unique_concepts:
        ##if a subsample of the concpets is used - comment in
        for concept in top_concepts:
            ##if the concept is in the dictionary, check the list of words
            if concept in con_dict:
                words = con_dict[concept]
                ##if there is only one word for the concept, append it to the language list
                if len(words)==1:
                    word_list_lang.append(words[0])
                ##if there are more words for the concept, join them to a string with - as eliminater
                else:
                    word_str = '-'.join(words)
                    word_list_lang.append(word_str)
            ##if the concept is not in the list, append 0 as placeholder
            else:
                word_list_lang.append("0")
        ##append everything to the dictionary
        nelex_dict[lang] = word_list_lang
    
    
    ##create sample dict for testing, get the first 40 keys. If you want to reduce the concept, do so in compute_distmtx/compute_dm
    #sample_dict = dict(it.islice(nelex_dict.items(),0,4))
    
    #return sample_dict
    
    return nelex_dict

def create_lodict():
    '''
    creates the lodict with the sounds and PMI scores
    :return lodcit: the dictionary with the sounds and the corresponding pmi scores
    '''
    
    ##reads the sounds from asjp present in nelex (39)
    f = open('input/sounds.txt')
    sounds = np.array([x.strip() for x in f.readlines()])
    f.close()

    ##reads the pmi matrix for nelex from tarakas method online pmi
    f = open('input/pmi_matrix_nelex.txt','r')
    #f = open('pmiTest.txt','r')
    l = f.readlines()
    f.close()
    logOdds = np.array([x.strip().split() for x in l],np.double)
    
    ##creates the lodict with the sounds and the pmi matrix
    lodict = dict()
    for i in xrange(len(sounds)):
        for j in xrange(len(sounds)):
            lodict[sounds[i],sounds[j]] = logOdds[i,j]
            
            
    return lodict  

def write_sim(f, overall_sim): 
    
    for langPair, simArray in overall_sim.items():
        #print langPair
        #print type(simArray)
        f.write(langPair[0]+","+langPair[1]+"\t")
        np.savetxt(f, simArray, delimiter=" ", newline="\t")
        f.write("\n")
    f.close()
    
    
if __name__ == '__main__':
    
    
    start = timeit.default_timer()
    
    ###########language tree##############
    ##read the data
    nelex_dict = read_nelex_original()
    ##get the log odds
    lodict = create_lodict()
    ##set the path to save the distance matrix
    out_path_mtx = "output/nelex/nelexDstMtx.phy"
    out_sim = "output/nelex/simMtx.phy"
    ##compute the distance matrix
    dist_dict, overall_simMtr = compute_distmtx.compute_dm(nelex_dict, lodict)
    fout = codecs.open(out_sim, "wb", "utf-8")
    ##write the similarity matris for the bootstrap
    write_sim(fout, overall_simMtr)
    #open the file and write the distance matrix
    fout1 = codecs.open(out_path_mtx,"wb","utf-8")
    compute_distmtx.writePhy(dist_dict, fout1)
      
    #########create the trees for the overall language tree################
    ##set the path to save the two trees
    out_nj = "output/nelex/NELexPMITree+nj.nwk"
    out_fastme = "output/nelex/NELexPMITree+fastme.nwk"
    ##reconstruct the trees with R and save them to a file
    reconstruct_trees.reconstruct_langTree(out_path_mtx, out_nj, out_fastme)
      
      
    ##################do bootstrap#########################
    ##do the traditional bootstrap with the similarity matrix, we don't have to compute the similarities again, which saves a lot of run time
    ##the matrices are directly written in a folder
    #traditional_bootstrap.create_bootstrapMtx(nelex_dict, overall_simMtr)
    traditional_bootstrap.create_bootstrapMtx(overall_simMtr)
    
    ####################create trees for the bootstrap replicates and the consensus tree########################
    in_path_mtx = "output/bootstrap/phylip/*.phy"
    ##output folders for the trees
    out_nj = "output/bootstrap/allBootstrapTrees+nj.nwk"
    out_fastme = "output/bootstrap/allBootstrapTrees+fastme.nwk"
    reconstruct_trees.reconstruct_bootstrap_trees(in_path_mtx, out_nj, out_fastme)
    ##create majority-rule consensus tree with PAUP*
    ##NJ
    out_consensus_nj = "output/bootstrap/consensusTree+nj.con.tre"
    create_consensus_paup.reconstruct_consensus(out_nj, out_consensus_nj)
    ##FastMe
    out_consensus_fastme = "output/bootstrap/consensusTree+fastme.con.tre"
    create_consensus_paup.reconstruct_consensus(out_fastme,out_consensus_fastme)



    
    stop = timeit.default_timer()
    print "end time",stop-start
    