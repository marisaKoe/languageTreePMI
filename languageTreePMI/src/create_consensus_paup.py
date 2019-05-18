'''
Created on 18.12.2017

@author: marisakoe
'''

import os
from subprocess import Popen

def reconstruct_consensus(treeFile, contreeFile):
    '''
    reconstruct a consensus tree with Paup for one single concept.
    :param treeFolder:the path to the tree folder
    :param concept:the name of the concept
    '''
    ##creates a temporary directory in this folder
    wdir = "output/temp/"
    
    ##returning a string representing the current working directory
    treePath = os.getcwd()+"/"+treeFile
    #print treePath
    #conFile = treeFile.split("/")
    #contree = conFile[0]+"/"+conFile[1]+"/consensus+nj.con.tre"
    
    print contreeFile
    ##create the skript for paup and save the trees in a file
    paupSkript = """#Nexus
    Begin paup;
    set incr=auto;
    gettrees file="""+treePath+""";
    contree /strict=no majrule=yes grpfreq=no showtree=no treefile="""+os.getcwd()+'/'+contreeFile+""" replace=yes;
    q;
    end;
    """
    #print paupSkript
    ##open the paup script
    with open(wdir+'/paupCommands.nex','w') as f:
        ##write
        f.write(paupSkript)
    ##open a terminal calling paup with the script as a command
    p = Popen("paup -n paupCommands.nex > /dev/null",shell=True,cwd=wdir)
    os.waitpid(p.pid,0)
    #if everything is done delete folder
    for f1 in os.listdir("output/temp"):
        os.remove(os.path.join("output/temp",f1))

if __name__ == '__main__':
    pass