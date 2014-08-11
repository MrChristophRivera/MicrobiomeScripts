#### Thsese are wrapper scripts for calling Quiime functions (I would rather do it this way)  
#This is kinda silly, in reality I would just rather call the code directly rather than having the shell do it. 

#import the needed modules
import sys
sys.path.append('/Users/christopherrichardrivera/Documents/PythonScripts')   #This appends stuff to the path this way I can access the function that I have written. 

import os
import subprocess as sp
from os.path import join

def assign_taxonomy(path, otus='otuFinal.fasta', outputdir = ""):
    #Uses Qiime to assign taxonomies. 
    #Inputs: path to the files, otus is the otu files
    #It calls suprocess to call rdp. 
    
    os.chdir(path)    #change the directory to the one of interest
    cmd = 'assign_taxonomy.py -i ' + join(path,otus) + ' -m rdp ' 
    if outputdir != "":
        cmd = cmd + '-o '+ join(path,outputdir)
    sp.call(cmd.split())
    print cmd

def assign_taxonomyParrallel(path, otus = 'otuFinal.fasta' , outpath = 'rdp_assigned_taxonomy', jobs = 2):
    #Usses the parrallel_assign_taxonomy_rdp to use the multicore  processors
    #Inputs 
    #   path: The path to all the files. 
    #   otu : the File name contining the outs (actually this could be any actual fasta file. 
    #   outpath: the name for the output data. Remember that the cmds will all be relative.
    #   jobs: The number of jobs to be run ( I believe that this should be less than the number of cores. Its unclear to me how much memory one needs. 
    
    cmd = 'parallel_assign_taxonomy_rdp.py  -i ' + join(path, otus) + ' -o ' + join(path,outpath) + ' -O ' + str(jobs)   
    sp.call(cmd.split())
    
def createOtuTable(path, OTUmap = 'seqs_otus.txt', output= 'otu_table.biom', taxonomy = 'rdp_assigned_taxonomy/otuFinal_tax_assignments.txt'  ):
    #Wrapper script to call the Qiime command for creating and otuTable and printing out summary information. It also uses the print biom_table summary to print out a summary.
    #Inputs 
    #   path: The path to all the files. 
    #   OTUmap: The name of the file that holds the otu mappings
    #   output: The name of the output file (This is a biom json format). 
    #   taxonomy : The location of the taxonomy file 
    
    #Make the Biom table
    cmd = 'make_otu_table.py -i ' + join(path, OTUmap) + ' -t ' + join(path, taxonomy) + ' -o ' + join(path,output)
    print cmd
    sp.call(cmd.split())
    

def printBiomSummary(path, biomTable = 'otu_table.biom', WRITE = False):       
    #print out the biom information 
    cmd = 'print_biom_table_summary.py -i ' +join(path,biomTable) 
    #print cmd
    #print sp.check_output(cmd.split())
    proc =sp.Popen(cmd.split(), stdout = sp.PIPE).communicate()[0]
    print proc
    
    #write it to file if true
    if WRITE == True:
        fout = open(join(path,biomTable.split('.')[0]) + '_results.txt', 'w')
        fout.write(proc)
        fout.close()
        
        
    
    return proc 
    
def assignTaxonomyAndCreateOtuTable(path, otus='otuFinal.fasta', OTUmap = 'seqs_otus.txt', output= 'otu_table.biom', taxonDir ='rdp_assigned_taxonomy'):
    """Uses qiime with rdp to assign taxonomies and create an otu Table"""
    #inputs:
    #   path: the path to the input files
    #   otus : Name of the fasta file that contains the outus to generate the otutable and the otu map. 
    #   OtuMap: Name of the tab delimited file that assigns reads to otus. Used to create an otu table. 
    #outputs:
    #   taxonDir: Name for the directory that holds the taxonmic information. Two files are stored here. 
    #   taxonomy: taxonomy is the name of the taxonomic text file from rdp. 
    
    assign_taxonomy(path, otus, taxonDir)
    
    #Create the taxonomy name
    taxonomy =  otus.split('.')[0] + '_tax_assignments.txt'
    taxonomy = join(taxonDir, taxonomy)    #name for the taxonomy text file to be fed into the createOtuTable function. 
    createOtuTable(path, OTUmap, output, taxonomy)
    Write= True
    printBiomSummary(path,output,Write)

###########################################################


def generateAlignmentAndTree(path, otus ='otuFinal.fasta', alignmentMethod ='muscle', treeMethod ='fasttree'):
    #Wrapper script to call Qiime to call several programs to generate a multiple sequence alignment and and a tree. 
    #Inputs:
    #    path: The root path wit the data
    #    otus: The fasta file that holds the rerpresentative otu files
    #    alignmentMethod: The program used for creating an alignment: 
    #        Possible methods are : '    pynast, infernal, clustalw, muscle, infernal, mafft 
    #    treeMethod: The program used to create a tree:
    #        Possible Methods are : 'clearcut, clustalw, fasttree_v1, fasttree, raxml_v730, muscle'
    # Note to self: Because I know more about phylogenetic trees I could use alternate methods. I want to get a computer for this. 
    
    #make a directory to place the results in.
    MSA = join(path, 'MSA_Tree')
    try: 
        os.makedir(join(MSA))
    except: 
        OSError
    
    #Generate the Alignment
    cmd =  'align_seqs.py -i ' + join(path, otus) + ' -m ' + alignmentMethod  + ' -o ' + MSA
    print cmd
    sp.call(cmd.split())
  
    #Generate the Tree   
    cmd2 = 'make_phylogeny.py -i ' +join(MSA,'otuFinal_aligned.fasta') + ' -t ' + treeMethod
    sp.call(cmd2.split())
    
    
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
    




if __name__== '__main__':
    ##path= '/Users/christopherrichardrivera/Documents/Pond_Project/Digester_Project/DigesterExperiment/Digester1'
    #path= '/Users/christopherrichardrivera/Desktop/Test'
    #assignTaxonomyAndCreateOtuTable(path, otus='otuFinal1.fasta', OTUmap = 'seqs_otus0_5.txt', output= 'otu_table1.biom', taxonDir ='rdp_assigned_taxonomy1')
    path = '/Users/christopherrichardrivera/Desktop/test'
    printBiomSummary(path, biomTable = 'otu_table4.biom', WRITE = True)
    
    
   