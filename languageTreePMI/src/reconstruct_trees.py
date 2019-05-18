'''
Created on 15.12.2017

@author: marisa
'''

import glob
import rpy2.robjects as r
from rpy2.robjects.packages import importr

#import the basis of R
base = importr("base")
utils = importr("utils")
stats = importr("stats")
#imports ape, required for phangorn
ape = importr("ape")
#imports phangorn
phangorn = importr("phangorn")

def reconstruct_langTree(in_path_mtx, out_nj, out_fastme):
    '''
    Reconstructs Neighbor Joining and FastMe language tree for the nelex
    The trees are reconstructed using R and the packages ape and phangorn
    '''

    t = utils.read_table(in_path_mtx, skip=1, row_names=1)
    mx = base.as_matrix(t)
    dm = stats.as_dist(mx)
    #nj trees
    tree = ape.nj(dm)
    ape.write_tree(tree, file=out_nj)
    #fastme trees
    tree1 = ape.fastme_bal(dm, nni=True, spr=True, tbr=False)
    ape.write_tree(tree1, file=out_fastme)
        
def reconstruct_bootstrap_trees(in_path_mtx, out_nj, out_fastme):
    '''
    reconstructs the bootstrap trees for all distance matrices in the folder
    the bootstrap trees are normal NJ and FastMe trees, computed of sampled data
    all trees are saved in one file, which can be read as input to PAUP* to construct a majoriy rule consensus tree
    :param in_path_mtx:the path for the distance matrices folder
    :param out_nj:the filename for all nj trees
    :param out_fastme:the filename for all fastme trees
    '''
    #print "in reconstruct trees"
    list_matrices = glob.glob(in_path_mtx)
    #print list_matrices
    for f in list_matrices:

        t = utils.read_table(f, skip=1, row_names=1)
        mx = base.as_matrix(t)
        dm = stats.as_dist(mx)
        #nj trees
        tree = ape.nj(dm)
        ape.write_tree(tree, out_nj, append=True)
        #fastme trees
        tree1 = ape.fastme_bal(dm, nni=True, spr=True, tbr=False)
        ape.write_tree(tree1, out_fastme,append=True)



if __name__ == '__main__':
    pass