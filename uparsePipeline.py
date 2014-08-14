# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 11:49:27 2014

@author: chrisrichardrivera
"""
#import required modules
import sys
import os
Path = '/Users/chrisrichardrivera/Documents/Code/Python/MicrobiomeScripts/'
sys.path.append(Path)  #This appends stuff to the path this way I can access the function that I have written. 
from UparseWrapperFunctions import *
from parseUC import *
import subprocess as sp
from os.path import join, split
from QiimeWrappers import *


#####The functions
def getNames(path):
    """returns a list of sample names from the files"""
    return [f.split('_')[0] for f in os.listdir(path) if f.endswith('.fastq') or f.endswith('.fasta')]    

    
def runUparse(path,Merge = True, truncqual = 5, minmergelen = 50, minovlen = 20, maxmergelen = 1000,maxee = 0.5, truncate=False, L=400, fastqout = False,otuid =0.97, ChimeraDatabase='/Users/chrisrichardrivera/Documents/BioinformaticDB/gold.fa', printCommand =False, description = ''):
    ''' Calls the Uparse Pipeline functions to run the Uparse Pipeline'''
    #Inputs:
    #      path: The path to the original fastq files that need to be merged'''
    #Merging    
    #      Merge: If true merge the files
    #      Merge: If true: merge the files
    #      truncqual: The qscore for whicht to truncate below- for merge
    #      minmergelen: minimum merge length - for merging
    #Filtering      
    #      minovlen:    minumum overlap length: for mergin
    #      maxee:    max expected errors for the given read, if read has greater than maxee discard. 
    #      truncate:    If true truncate the reads to given length L
    #      L:     lengh tpo truncate to
    #      fastqout: if True, output fastq files.           
    #Clustering
    #      otuid: The threshold for clustering. 
    #Chimera removal
    #      ChimeraDatabase: The absolute path the chimera database that one is using-the gold fa
    #All
    #       printCommnad tells the functions to print out the commands.  
    
    parent = split(path)[0]   #get the parent path
    
    #if Merge is true, merge the fastq present, and place the merged files in the new directory. 
    if Merge == True:
        
        mergePath = join(parent, 'Merge' + description)   #create the path to the merge directory
        mergeAllReads(path, mergePath ,truncqual, minmergelen, minovlen, maxmergelen, printCommand)
        fastqPath = mergePath
    else:
        fastqPath = path   #if done twant to merge, need to tell the program where to go. 
    
    #Quality filter the reads
    qualityPath = join(parent, 'Filtered' + description)   #the target for the next step as well.
    qualityFilterAllFastq(fastqPath, qualityPath, maxee, truncate, L, fastqout,printCommand)
    
    #make a new child directory that can hold the results for the remaining files. I will keep all the files. 
    results = join(parent, 'Results'+ description )
    
    try:
            os.mkdir(results)
    except:
            OSError

    #Rename the headers for the sequences and combine the files into one large file. 
    #Note to self, we have the results as FNA files. This makes it easier.. 
    names = getNames(qualityPath)   #get the names for the next step
    index = 1    #this is the starting index
    Combined = join(results,'Combined' + description + '.fna' )    #Set up the file name for the new file. It will be in results. 
    
    
    renameAndCombineFastas(qualityPath, names, index, Combined)   #rename and combine the files. 
    
    #Dereplicate the reads
    Dereplicated = join(results, 'Dereplicated' + description + '.fna')    #Set up the path to the dereplicated file. 
    dereplicateReads(Combined, Dereplicated,printCommand)
    
    #Abundance Sort the Reads
    Sorted = join(results, 'Sorted' + description + '.fna')              #set up the path to sorted target file. 
    minsize = 2
    abundanceSortReads(Dereplicated, Sorted,minsize, printCommand)
    
    #Cluster the reads
    Clustered = join(results, 'Clustered' + description + '.fna')       #Set up the path to the clustered file. 
    clusterOtus(Sorted, Clustered, otuid,printCommand)             #Cluster the reads with an otuid threshold #Keep greater than 0.97 if you can. 
    
    #Remove remaining Chimeras
    Otus = join(results, 'Otus' + description + '.fna')
    removeChimeras(Clustered, Otus, ChimeraDatabase,printCommand)
    
    #Rename the reads
    renameReads(Otus, 'denovo') 
    
    #Map the sequences back to an UC file for down stream analysis. 
    OtuUC = join(results,'Otus' + description + '.uc')
    mapSeqToOTUS(Combined, Otus, OtuUC, otuid,printCommand)   #map the sequences in the combined file ont to the Otus.fna file to create an outuUC file with threshold otuId. 
    
    #Parse the UC file
    OtuTxt = join(results, 'Otus' + description + '.txt' + description)    #Path to the text file. 
    parseUCandWriteOtus(OtuUC, OtuTxt)
    

if __name__ == '__main__':
    
    #Set up the parameters
    path= '/Users/chrisrichardrivera/Desktop/test/sequences'
    Merge = True
    truncqual = 5
    minmergelen = 50
    minovlen = 20
    maxmergelen = 1000
    maxee = 0.5
    truncate=False
    L=400
    fastqout = False
    otuid =0.97
    ChimeraDatabase='/Users/chrisrichardrivera/Documents/BioinformaticDB/gold.fa'
    printCommand = True
    description = 'Test'
    runUparse(path,Merge, truncqual, minmergelen, minovlen, maxmergelen, maxee, truncate, L, fastqout ,otuid, ChimeraDatabase, printCommand, description)
#    fin = '/Users/chrisrichardrivera/Desktop/ResultsTest/OtusTest.uc'
#    fout = '/Users/chrisrichardrivera/Desktop/ResultsTest/OtusTest.txt'
#    
#    parseUCandWriteOtus(fin ,fout)
#    