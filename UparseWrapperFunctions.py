# Load the desired modules
import sys

#import parseUC
import os
import subprocess as sp
from os.path import join
import gzip
#from time import time

sys.path.append('/Users/chrisrichardrivera/Documents/PythonScripts/MicrobiomeScripts')   #append the path to libraries. 
from QiimeWrappers import *


##################################
#Merging

def mergeReads(forward, reverse, outputName ='merged.fastq' ,truncqual = 5, minmergelen = 50, minovlen = 20, maxmergelen = 1000, printCommand=False):
    #pass
    '''Call the usearch fastq_mergepaired program'''
    #Inputs :    
    #            Forward: path to the fastq file with the forward reads. 
    #            Reverse: path to the fastq file with the reverse reads.
    #            truncqual: qvalue below which the sequences is truncated
    #            minovlen: minimum overlap for the merged pairs. If below threshold, read pair discarded
    #            minmergelen: minimum length of the paired read after merging. If below the minimum length, read pair discarded. 
    #            maxmergelen: maxiumum length of the paired read after merging. If above the maximum Length, read pair discarded
    #            printCommand: If True, print out the commmand sent to sp.call for calling usearch. 
    
    #Set up the commands and Use subprocess to call usearch.
    cmd = 'usearch -fastq_mergepairs ' +  forward + ' -reverse ' + reverse + ' -fastq_truncqual  ' + str(truncqual) + ' -fastq_minmergelen ' + str(minmergelen) + ' -fastq_maxmergelen ' + str( maxmergelen) + ' -fastq_minovlen ' + str(minovlen) + ' -fastqout ' +  outputName                                   
    
    if printCommand == True:
        print cmd
    sp.call(cmd.split())
   
   
def mergeAllReads(path, outputPath ,truncqual = 5, minmergelen = 50, minovlen = 20, maxmergelen = 1000, printCommand=False):
    '''Marges all fastq files in a given directory'''
    # Inputs:
    #       path: Path to the input files
    #       outputPath : path to the output path. 
    
    #Identify all the foward and reverse fastq    
    forward = [f for f in os.listdir(path) if f.endswith('R1_001.fastq')]
    reverse = [f for f in os.listdir(path) if f.endswith('R2_001.fastq')]
           
    #If the directory for the output files does not exist, make it. 
    try:  
        os.mkdir(outputPath)
    except:
        OSError
    
    #Sequentially merge the files. 
    for i in range(len(forward)):
        #get the path to the file. 
        f = join(path,forward[i] )
        r = join(path, reverse[i] )
        prefix = forward[i].split('.')[0]   #get the prefix
        
        #Make the name
        Merge = join(outputPath,prefix + '.fastq' )
        #merge the reads
        
        mergeReads(f, r, Merge,truncqual, minmergelen, minovlen,maxmergelen, printCommand)

#########################################################
#Qualiy Filtering. 


def qualityFilterFastq( fastq, fasta='filtered.fasta', maxee = 0.5, truncate=False, L=400, fastqout = False, fastqoutput ='filtered.fastq', printCommand = False):
    """Calls usearch -fastq_filter to qaulity fitler reads in a fastq file"""
       
    #Inputs:     
    #   fastq: Absolute path to the fastq file.        
    #   fasta: Absolute path to the fastq file 
    #   maxee: Maximum expected errors, if maxee is greater than threshold, discard read. 
    #   truncate: If true, truncate the read to given length L. 
    #   L: Truncation length. Truncates to Length L, if less than L, discard read. 
    #   fastqout: If true, create a fastq file as well.
    #   fastqputput: the path to the fastq output
    #   printCommand: If true print the command. 
    
    #set up the base command. 
    cmd = 'usearch -fastq_filter ' + fastq + ' -fastq_maxee ' + str(maxee) + ' -fastaout ' + fasta
    
    #add truncation command if needed.     
    if truncate==True:
        cmd += ' -fastq_trunclen ' + str(L)   #The base command
    
    #if the a fastq is desired just add the addiional command to the sting 
    if fastqout == True:
        cmd +=' -fastqout ' + fastqoutput    #The command if a fastq is desired. 
        
    #use the subprocess program to use the shell to call the uparse function.     
    sp.call(cmd.split())
    
    #print the cmd if desired
    if printCommand ==True:
        print cmd

def qualityFilterAllFastq(path, outputpath, maxee = 0.5, truncate=False, L=400, fastqout = False, printCommand = False):
    '''Calls usearch to quality filter all fastq files within the target path. '''

    #Inputs:
    #    Path: path to the Fastq files
    #    outputpath: path for the output fasta files.
    #    See qualityFilterFastq for more details on the inputs. 
    
    fastqs = [f for f in os.listdir(path) if f.endswith('fastq')]    #get the fastq file names. 
    fastas = [join(outputpath, f.split('.')[0] + '.fasta') for f in fastqs]    #set the output names for the fastqs
    
    # generate names for the fastqoutput files, may or may not be needed. 
    fastqOutputPath = outputpath + 'Fastq'
    fastqouts = [join(fastqOutputPath,f) for f in fastqs]
    
    #generate the path for the input fastqfiles
    fastqs = [join(path, f) for f in fastqs]
    
    
    #Make directories for output if desired. 
    try:
        os.mkdir(outputpath)
    except:
        OSError
    
    #Make output directory for fastqfiles if needed and the fastqouts if neeeded. 
    if fastqout == True:
        try:
            os.mkdir(fastqOutputPath)
        except:
            OSError
 
    #call the usearch commmand for each
    [qualityFilterFastq(fastqs[i],fastas[i], maxee, truncate, L,fastqout,fastqouts[i], printCommand ) for i in range(len(fastqs))]
    


############################################
#Rename the reads and fasta files. 

                                                                                                                                                          
def renameReads(fasta, sampleName='Sample', index = 1):
    #This Reads of a fasta file. It renames the reads as >Sample_index eg. >Sample_1, >Sample_2
    #Inputs:
    #       fasta = the name of the fasta file
    #       sampleName = the name of the sample prefix
    #       index =   the starting index.
    lines = [line for line in open(fasta)]   #get the name
    
    
    #Rename the Header Lines
    for i in range(len(lines)):
        if lines[i].startswith('>'):
            lines[i] = '>' + sampleName + '_' + str(index) + '\n'
            index +=1
    f_out = open(fasta, 'w')
    f_out.writelines(lines)
    f_out.close()
    
def combineFastas(path, outputName = 'combined.fasta'):
    """Comines the fasta files in a directory into to one file"""
    #Inputs:    
    #    path: The directory with the fasta files
    #    outputName : The path to the fasta file. 
  
    files = [join(path, f) for f in os.listdir(path) if f.endswith('fasta')]    #get a list of the files. 
     #combined = open(outputName, 'w')
    
    lines = []
    for f in files:
        lines += open(f).readlines() 
    #open a new file and write the fasta
    fout = open(outputName, 'w')
    fout.writelines(lines)
    fout.close()
    
  
def renameAndCombineFastas(path, names, index = 1, output='combined.fasta'):
    """Ranmes and combines fastq names""" 
    #Inputs:
    #   path: The path to the files
    #   names: The sample name to be applied to each sequence
    #   index: The starting number for the indexes to be passed to renameReads
    #   outputName: Obvious
    
   
    files = [join(path,f) for f in os.listdir(path) if f.endswith('fasta')] 
    
    #error check
    if len(names) != len(files):
        print 'The number of label names does not equal the number of Fastas'
    else:
        #rename the files
        [renameReads(files[i], names[i], index) for i in range(len(files))]
        #combine the files: 
        combineFastas(path, output)


############
# Dereplicate Reads, abundnace sort and cluster reads.
def dereplicateReads(fasta, output ='dereplicated.fasta', printCommand=False):
    #Wrapper function for Usearch that dereplicates reads. 
    #Inputs:
    #   fasta: the target merged, filtered, renamed and combined fasta
    #   output: a file name 
    #   printCommand: if true print Command
    
    cmd = 'usearch -derep_fulllength ' + fasta + ' -output ' + output + ' -sizeout'
    #print cmd.split()
    sp.call(cmd.split())
    
    if printCommand == True:
        print cmd
        
def abundanceSortReads(fasta, output = 'abundanceSorted.fasta', minSize = 2, printCommand=False):
    #Wrapper function for Usearch that abundance sorts Reads and removes clusters that are less than a certain size. 
    #Inputs: 
    #   fasta: the target fasta to be dereplicated
    #   output: the output file name
    #   minSize: The minimum size of the clusters
 
    cmd = 'usearch -sortbysize ' + fasta + ' -output ' + output + ' -minsize ' + str(minSize)
    sp.call(cmd.split()) 
    if printCommand == True:
        print cmd
  
def clusterOtus(fasta, output ='otus.fasta', otuid =0.97,printCommand=False):
    #Wrapper function for Usearch clusters otus and and removes chimeras.
    #Note to self, this function bassicaly removes chimeras, and then picks a rep for a cluster. It does not have any real size info that is important. 
    #Inputs: 
    #   fasta: the target fasta to be dereplicated
    #   output: the output file name
    #   otuid:sets the identity between the cluster and the sequences. Default is 0.97, it is recomended that 0.97 is minimum.
    
    cmd = 'usearch -cluster_otus ' +  fasta + ' -id ' + str(otuid)  + ' -otus ' + output
    sp.call(cmd.split())
    if printCommand == True:
        print cmd
    
def removeChimeras(fasta, output, database, printCommand=False):
    #Wrapper function for Usearch that does one additional step of chimera removal
    #This function uses usearch align the otu sequencs against the reference database and determine if they are chimeras. 
    
    #Inputs: 
    #   fasta: the fasta file with the reference otu sequences. 
    #   output: the output file name
    #   database: path to the reference database for which the sequences are compared. 
    
    cmd = 'usearch -uchime_ref ' + fasta + ' -db ' + database + ' -nonchimeras '  + output + ' -strand plus' 
    sp.call(cmd.split())
    if printCommand == True:
        print cmd

def mapSeqToOTUS(fasta, otuFasta, output= 'otus.uc', identity= 0.97 , printCommand=False):
    #Wrapper function for Usearch that uses the global usearch option to cluster the sequences in a file and the cluster them against reference otus. 
    #This outputs a uc file, which can later be parsed to make a text file. 
    #This function uses usearch align the otu sequencs against the reference database and determine if they are chimeras. 
    
    #Inputs: 
    #   fasta: the fasta file with  original sequences. 
    #   outFasta is the fasta file that contains the otu representatives. 
    #   output: the output file name, The output file is in uc format- its a tab delimited file 
    #   database: path to the otu sequences to which the sequences are compared.
    #   id: the %id threshold 
    
    cmd ='usearch -usearch_global ' + fasta + ' -db ' + otuFasta +  ' -strand plus ' + ' -id ' +str(identity) + ' -uc ' +  output
    sp.call(cmd.split())
    if printCommand == True:
        print cmd





#######Misc functions. 

def cropPairedReads(read1,read2, inputpath, outputpath, pairedpath, HeadCrop = 30, MinLen=200):
    """ Calls trimmomatic to head trim the reads for simple processing."""
    #note this makes 4 different files for the paried reads. This is really the best way to do it.. 
    #read1, read2 are the names of the paired reads, inputpath is the name of the directory, output path is the name for the "unpaired reads" and paired path is the name for the "paired reads"
    #The paired reads are reads for which both member survived the cutting. 
    #Make an output directory
    try:
        os.mkdir(outputpath)
    except:
        OSError
        
    try:
        os.mkdir(pairedpath)
    except:
        OSError    
    
    cmd = 'java -jar /Users/christopherrichardrivera/Trimmomatic-0.32/trimmomatic-0.32.jar PE ' + join(inputpath,read1) + ' ' + join(inputpath,read2) +' ' + join(pairedpath, read1) + ' ' + join(outputpath,read1) + ' ' + join(pairedpath, read2) + ' ' + join(outputpath,read2) + ' ' +  'HEADCROP:'+ str(HeadCrop) + ' ' + 'MINLEN:' + str(MinLen)
    #return cmd
    sp.call(cmd.split())




def decompressGZ(f_in, savepath=''):
    """Calls gzip to open a file and decompress it. It renames the file."""
    #Note this assumes that the file is relative path not absolute.  
    
    f_out = f_in.replace('.gz', '')
    
    if savepath =='':    
        f_out = join(savepath, f_out)    #
        print f_out
    else:
        f_out = join(savepath, f_out.split('/')[-1])
    #Open the file with gzip, read and write it in noncompressed format.
    f_in = gzip.open(f_in)
    f_out = open(f_out, 'w')
    f_out.writelines(f_in)
    
    #close the unneed files
    f_in.close()
    f_out.close()
 
            
######################################################

def countSeqFast(fname):
    """ Uses the command line to count the number of lines"""
    cmd = ['wc', '-l' ,fname]
    proc = sp.Popen(cmd, stdout = sp.PIPE)
    size = proc.communicate()[0].split()
    
    if fname.endswith('fastq'):
        return int(size[0]) /4
    if fname.endswith('fasta') or fname.endswith('fa'):
        return  int(size[0])/2 
    return None
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
if __name__ == '__main__':
    
    #path= '/Users/christopherrichardrivera/Documents/Pond_Project/Digester_Project/DigesterExperiment/Digester4Trimmed'
    #tabulateReads(path)
    
    
    