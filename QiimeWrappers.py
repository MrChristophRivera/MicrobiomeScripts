# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 11:22:04 2014

@author: chrisrichardrivera
"""

#import the required modules

import os
import subprocess as sp
from os.path import join,split

#Wrapper functions
def assign_taxonomy(otus, method = 'rdp', printCommand=False):
    '''Calls Qiime to call RDP to assign Taxonomy'''
    #Inputs:
    #   otus: absoulte path tho the otu fasta file
    #The fuction places the output to the directory one level up. 
    parent = split(otus)[0]   #get the parent directory
    child = join(parent, 'rdp_assigned_taxonomy')
    
    cmd = 'assign_taxonomy.py -i ' + otus + ' -m ' + method + ' -o ' + child

    sp.call(cmd.split())
    if printCommand == True:
        print cmd
    
def createOtuTable(mapfile, taxonomy, printCommand=False):    
    #Calls Qiime python script to generate an OTU table
    #Inputs:
    #     mapfile:absolute path to the map text file that maps otus to sequences
    #     taxonmy: absolute path to the map text file that maps otus to sequences.  
    
    biom = mapfile.split('.')[0] +'.biom'    #set up a path for the biom. 
    #Make the Biom table
   
    cmd = 'make_otu_table.py -i ' + mapfile + ' -t ' + taxonomy + ' -o ' + biom
    sp.call(cmd.split())
    
    #if printCommand is True: print the command
    if printCommand == True:
        print cmd
    
def writeBiomSummary( biom = 'otu_table.biom'):       
    #writes a summary of a Biom file
    summary = join(split(biom)[0], 'biomSummary.txt' )  #set up a name for the summary
    cmd = 'biom summarize-table -i ' + biom + ' -o ' + summary  #set up the command
    sp.call(cmd.split())     #call the command

#function that calls the scripts.     
def assignTaxonomyAndCreateOtuTable(path, printCommand=True):
    """Uses qiime with rdp to assign taxonomies and create an otu Table"""
    #inputs:
    #   otus: The path to the directory that holds the otu file of interest.  
    #   printCommand: If true print the commands
    otus = [join(path,f) for f in os.listdir(path) if f.startswith('Otus') and f.endswith('.fna')][0] #get the path to the fasta file with the sequences 
    method = 'rdp'

    assign_taxonomy(otus, method, printCommand)

    
    #create an otu table
    mapfile = [join(path,f) for f in os.listdir(path) if f.startswith('Otus') and f.endswith('.txt')][0]  #get the absolute path to the map file. 
    rdpPath = join(path, 'rdp_assigned_taxonomy')   #string for path to the rdp taxonomy file
    taxonomy = [join(rdpPath, f) for f in os.listdir(join(rdpPath )) if f.endswith('tax_assignments.txt')][0]   #get the rdp taxonompy files
    createOtuTable(mapfile, taxonomy, printCommand=False)
    
    #write a Biom summary table
    biom = mapfile.split('.')[0] + '.biom'   #get the name for the biom file
    writeBiomSummary(biom)

###########################################################


def generateAlignmentAndTree(otus, alignmentMethod ='muscle', treeMethod ='fasttree', printCommand = False):
    #Wrapper script to call Qiime to call several programs to generate a multiple sequence alignment and and a tree. 
    #Inputs:
    #    otus: Tabsoulute path to the otus file. 
    #    alignmentMethod: The program used for creating an alignment: 
    #        Possible methods are : '    pynast, infernal, clustalw, muscle, infernal, mafft 
    #    treeMethod: The program used to create a tree:
    #        Possible Methods are : 'clearcut, clustalw, fasttree_v1, fasttree, raxml_v730, muscle'

    
    #make a directory to place the results in.
    parent = split(otus)[0]
    MSA = join(parent, 'MSA_Tree')
    
    try:  
        os.makedir(MSA)
    except: 
        OSError
    
    #Generate the Alignment
    cmd =  'align_seqs.py -i ' + otus + ' -m ' + alignmentMethod  + ' -o ' + MSA
    print cmd
    sp.call(cmd.split())
  
    #Generate the Tree 
    aligned = [join(MSA,f) for f in os.listdir(MSA) if f.endswith('.fasta')][0]
    cmd2 = 'make_phylogeny.py -i ' +aligned + ' -t ' + treeMethod
    sp.call(cmd2.split())
    if printCommand == True:
        print cmd
        print cmd2


#rarefaction analysis.     
def performRarefactionAndAlphaDiverstiy(path, Min = 10 , Max = 110, Step = 20 ):
    #Wrapper script that uses Qiime scripts to perform rarefaction and conduct the most basic alpha diversity analysis 
    #Inputs:
    #   path: The directory holding the biom file. We are expecting that I have not changed these names
    #   Min : the minimum that everything will be rarefied to 
    #   Max: The maximum that every sample will be rarefied to
    #   Step: The step size
    
    #Note its very important to always use absolute paths. 
    
    
    #Rarefy the Data
    cmd = 'multiple_rarefactions.py -i ' + join(path, 'otu_table.biom') + ' -m ' + str(Min) + ' -x ' + str(Max) + ' -s ' + str(Step) + ' -o '  + join(path,  'rarefied_biom/')
    sp.call(cmd.split())
    
    #Calculate alpha diversities
    cmd = 'alpha_diversity.py -i ' + join(path, 'rarefied_biom/') + ' -t ' + join(path, 'MSA_Tree', 'otuFinal_aligned.tre') + ' -o ' + join(path, 'alpha_diversity/')
    sp.call(cmd.split())
    
    #Collate the alpha diversity information
    cmd = 'collate_alpha.py -i ' + join(path, 'alpha_diversity/') + ' -o ' + join(path,'alpha_summary/')
    sp.call(cmd.split())
    
######Depricated
def assign_taxonomyParrallel(path, otus = 'otuFinal.fasta' , outpath = 'rdp_assigned_taxonomy', jobs = 2):
    #Usses the parrallel_assign_taxonomy_rdp to use the multicore  processors
    #Inputs 
    #   path: The path to all the files. 
    #   otu : the File name contining the outs (actually this could be any actual fasta file. 
    #   outpath: the name for the output data. Remember that the cmds will all be relative.
    #   jobs: The number of jobs to be run ( I believe that this should be less than the number of cores. Its unclear to me how much memory one needs. 
    
    cmd = 'parallel_assign_taxonomy_rdp.py  -i ' + join(path, otus) + ' -o ' + join(path,outpath) + ' -O ' + str(jobs)   
    sp.call(cmd.split())



if __name__== '__main__':
    
    
    fin= '/Users/chrisrichardrivera/Desktop/test/ResultsTest/OtusTest.txt'
    assign_taxonomy(fin)