#### These are master scripts that will use some of the scripts from modules that I have written to more automate the pipeline. 

#import the needed modules
import sys
sys.path.append('/Users/chrisrichardrivera/Documents/PythonScripts/MicrobiomeScripts')  #This appends stuff to the path this way I can access the function that I have written. 
from UparseWrapperFunctions import *
from parseUC import *
import os
import subprocess as sp
from os.path import join
import gzip
from QiimeWrappers import *




def runUparseToPickOtusAndMapReads( path, sampleNames,  gz =True, forwardID= 'R1', reverseID= 'R2', minmergelen = 50, truncqual = 5, maxee = 0.5, refdb ='/Users/christopherrichardrivera/Documents/BioinformaticDB/gold.fa'  ,ucName = 'otu.uc', seq_otuName = 'seqs_otus.txt'): 
    #This client funciton runs the entire UparsePipeline
    #See the code for the api functions that this client calls to understand all the parameters. 
    #This function calls multipe functions to use UParse to pick otu representative otu sequences for reasds.
    
    
    SequenceCount = []    #This is to store all the sequence counts. So That we Know how the Pipeline Effects it 
    #Step = ['Total Reads', 'Merged Reads', 'Filtered Reads', 'Dereplicated Reads', 'Sorted Reads', 'Clustered Otus', "Otus with Chimeras Removed"]           #This is to store strings so I can have a simple annotation for Printing later. 

    #Set up that names and paths for all the files to be used in the downstream functions.  
    outputPath = 'pairedFiltered'
    tempPath ='tempFasta'    
    Path = join(path, outputPath)           #set up the path to the stuff
    SeqsFasta = join(path, 'reads.fasta')    #This is the file that contains all the merged and filtered files.
    otus = join(path,'otus.fasta')   
    otuNoChimera = join(path, 'otuFinal.fasta')
    uc = join(path, ucName)   
    seq_otus = join(path, seq_otuName)       #The tab delimited file that is used by Qiime to make an otu file.     
                    
    #Merge and Filter Paried End Reads 
    (totalCount, mergeCount, filterCount) = mergeAndQualityFilterFastQFiles(path,gz, forwardID, reverseID, outputPath, tempPath, truncqual , minmergelen, maxee)
    ######return sequenceCount
    
    #Update the count and sequences
    SequenceCount.append(totalCount)
    SequenceCount.append(mergeCount)
    SequenceCount.append(filterCount)
    #Rename the sequences headers of the paired reads based on a sample name and combine them into one fasta. 
    
    index = 1
    renameAndCombineFastas(Path, sampleNames, index, SeqsFasta)
    
    #Dereplicate the reads
    dereplicateReads(SeqsFasta, otus)
    SequenceCount.append(countSeq(otus))
    
    #Abundance sort the reads and discard cluster less than 2  
    abundanceSortReads(otus,otus)
    SequenceCount.append(countSeq(otus))
    
    #OTU cluster the reads and remove chimeras
    clusterOtus(otus, otus)
    SequenceCount.append(countSeq(otus))
    
    #Perform seconday Chimera Removal against a reference database. 
    removeChimeras(otus, otuNoChimera, refdb)
    SequenceCount.append(countSeq(otuNoChimera))
    
    #rename the otu sequences. 
    renameReads(otuNoChimera, 'denovo')
    
    #Map the sequences to the Otus and to make an otuFile. 
    mapSeqToOTUS(SeqsFasta, otuNoChimera ,  uc, identity = 0.97)
    
    #parse the UC file to make a tab delimited file that can be used to map the files
    parseUCandWriteOtus(uc, seq_otus)
    
    return SequenceCount

def printCountStats(SequenceCount):
    #This is a for loop that just prints the count stats so I can look at it in its glory
    Step = ['Total:  ', 'Merged:', 'Filtered:', 'Dereplicated:', 'Sorted:', 'Clustered:', "Chimeras Removed:"]  
    
    print 'Step'.ljust(20) + '\t' + 'Read Count' + '\t' +'Fraction'
    for i in range(len(SequenceCount)):
        print Step[i].ljust(20) + '\t' + str(SequenceCount[i]).ljust(10)  + '\t' +  str(SequenceCount[i]/float(SequenceCount[0]) *100)

def chartCountStats(SequenceCount):
    #This function will just plot a simple histogram showing how the number of reads changes with each step in creating otus. 
    pass
    
def getFastqNames(path):
    """returns a list of sample names from the files"""
    return [f.split('_')[0] for f in os.listdir(path) if f.endswith('.fastq')]    
    



################################################################################################################################################################################################ variants

def uparsePostMerge( path, sampleNames, maxee = 0.5, refdb ='/Users/christopherrichardrivera/Documents/BioinformaticDB/gold.fa'  ,ucName = 'otu.uc', seq_otuName = 'seqs_otus.txt'): 
    """ Simpler variant of the pipeline that does not do the merging (which has presumably already been done"""
    
    #Set up that names and paths for all the files to be used in the downstream functions.  
    outputPath = 'pairedFiltered'
    tempPath ='tempFasta'    
    Path = join(path, outputPath)           #set up the path to the stuff
    SeqsFasta = join(path, 'reads.fasta')    #This is the file that contains all the merged and filtered files.
    otus = join(path,'otus.fasta')   
    otuNoChimera = join(path, 'otuFinal.fasta')
    uc = join(path, ucName)   
    seq_otus = join(path, seq_otuName)       #The tab delimited file that is used by Qiime to make an otu file.     
                    
    #Filter Paried End Reads 
    fillterFastQFiles(path, outputPath,tempPath,  maxee)
    
    index = 1
    renameAndCombineFastas(Path, sampleNames, index, SeqsFasta)
    
    #Dereplicate the reads
    dereplicateReads(SeqsFasta, otus)
    
    #Abundance sort the reads and discard cluster less than 2  
    abundanceSortReads(otus,otus)
    
    #OTU cluster the reads and remove chimeras
    clusterOtus(otus, otus)
    
    #Perform seconday Chimera Removal against a reference database. 
    removeChimeras(otus, otuNoChimera, refdb)
    
    #rename the otu sequences. 
    renameReads(otuNoChimera, 'denovo')
    
    #Map the sequences to the Otus and to make an otuFile. 
    mapSeqToOTUS(SeqsFasta, otuNoChimera ,  uc, identity = 0.97)
    
    #parse the UC file to make a tab delimited file that can be used to map the files
    parseUCandWriteOtus(uc, seq_otus)
    
###################################################### (This lacks count of sequences)

def runUparsePipelineNoCounting( path, sampleNames,  gz =False, minmergelen = 50, truncqual = 5, minovlen = 10, maxee = 0.5, refdb ='/Users/christopherrichardrivera/Documents/BioinformaticDB/gold.fa'  ,ucName = 'otu.uc', seq_otuName = 'seqs_otus.txt'): 
    #This client funciton runs the entire UparsePipeline
    #See the code for the api functions that this client calls to understand all the parameters. 
    #This function calls multipe functions to use UParse to pick otu representative otu sequences for reasds.
    #its is variant that does does no counting. TO speed up the analysis. 
    #note that this has one more filtering parameter.. 

    #Set up that names and paths for all the files to be used in the downstream functions.  
    outputPath = 'pairedFiltered'
    tempPath ='tempFasta'    
    Path = join(path, outputPath)           #set up the path to the stuff
    SeqsFasta = join(path, 'reads.fasta')    #This is the file that contains all the merged and filtered files.
    otus = join(path,'otus.fasta')   
    otuNoChimera = join(path, 'otuFinal.fasta')
    uc = join(path, ucName)   
    seq_otus = join(path, seq_otuName)       #The tab delimited file that is used by Qiime to make an otu file.     
    forwardID= 'R1'         
    reverseID= 'R2'       
    #Merge and Filter Paried End Reads 
    #mergeAndQualityFilterFastQFiles(path,gz, forwardID, reverseID, outputPath, tempPath, truncqual , minmergelen, maxee) Old commmand
    mergeAndQualityFilterFastQFiles(path, gz, forwardID, reverseID , outputPath, tempPath,  truncqual, minmergelen, maxee, minovlen)
            
    index = 1
    renameAndCombineFastas(Path, sampleNames, index, SeqsFasta)
    
    #Dereplicate the reads
    dereplicateReads(SeqsFasta, otus)
    
    #Abundance sort the reads and discard cluster less than 2  
    abundanceSortReads(otus,otus)
    
    #OTU cluster the reads and remove chimeras
    clusterOtus(otus, otus)
    
    #Perform seconday Chimera Removal against a reference database. 
    removeChimeras(otus, otuNoChimera, refdb)
    
    #rename the otu sequences. 
    renameReads(otuNoChimera, 'denovo')
    
    #Map the sequences to the Otus and to make an otuFile. 
    mapSeqToOTUS(SeqsFasta, otuNoChimera ,  uc, identity = 0.97)
    
    #parse the UC file to make a tab delimited file that can be used to map the files
    parseUCandWriteOtus(uc, seq_otus)
    



                    
if __name__ == '__main__':
    gz = False   #tells the pipeline not to decompress
    minmergelen = 300
    truncqual = 5
    minovlen = 75
    maxee = 0.5
    refdb ='/Users/christopherrichardrivera/Documents/BioinformaticDB/gold.fa'
    ucName = 'otu.uc'
    seq_otuName = 'seqs_otus.txt'
    path = '/Users/christopherrichardrivera/Desktop'
    trimmed=['Test']

    for i in range(len(trimmed)):
        
        inputPath = join(path, trimmed[i])
        forward = [f for f in os.listdir(inputPath) if f.endswith('R1_001.fastq')]
        names = [f.split('_')[0] for f in forward]
        
    
        #assign the taxonomy and create a otu table
        assignTaxonomyAndCreateOtuTable(inputPath)
    
    