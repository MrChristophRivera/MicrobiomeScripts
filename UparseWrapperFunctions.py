# Load the desired modules
import sys
sys.path.append('/Users/chrisrichardrivera/Documents/PythonScripts/MicrobiomeScripts')   #This appends stuff to the path this way I can access the function that I have written. 

#import parseUC
import os
import subprocess as sp
from os.path import join
import gzip
#from time import time
from QiimeWrappers import *


#Below are my own wrapper scripts to be run so that I can use Uparse

def mergePairedEnds(forward, reverse, outputName ='merged.fastq' ,truncqual = 5, minmergelen = 50, minovlen = 20, maxmergelen = 1000):
    #pass
    """Wrapper function that calls the usearch fastq_mergepaired commands"""
    #Inputs :    forward is the forward sequence
    #            reverse is the reverse sequence
    #            truncqual is  the qvalue below which the sequences is truncated
    #            minovlen k is the minimum overlap for the merged pairs. The optimal value depned on the settings 
    #            minmergelen L is the minimum length of the paired read after merging. 
    #            maxmergelen L is the maxiumum length of the paired read after merging. 
    
    #Set up the commands and Use subprocess to call usearch.
    cmd = 'usearch -fastq_mergepairs ' +  forward + ' -reverse ' + reverse + ' -fastq_truncqual  ' + str(truncqual) + ' -fastq_minmergelen ' + str(minmergelen) + ' -fastq_maxmergelen ' + str( maxmergelen) + ' -fastq_minovlen ' + str(minovlen) + ' -fastqout ' +  outputName                                   
    #print cmd
    sp.call(cmd.split())
    #return cmd.split()
    
def qualityFilterReads( merged, outputName, maxee = 0.5):
    """Wrapper function that calls usearch to quality filter reads. If the read has an greater than 0.5 expected errors it is discareded."""       
    #Inputs:     merged is the input file name       
    #            maxee is the maximum expected errors (remember that the expected errors is simply the summation of all the probabilities of errors
    #            outputName is the name of the output file
    
    #Setup the usearch command and use subprocess to call it. 
    cmd = 'usearch -fastq_filter ' + merged + ' -fastq_maxee ' + str(maxee) + ' -fastaout ' + outputName     
    sp.call(cmd.split())


def qualityFilterAndTruncateFastq(fastq, outputName, maxee= 0.5, L= 400, fastqout = False, fastqOutputName=''):
    '''Wrapper function that calls usearch to quality filter reads and truncate the read length to some N length
    if the fastaout tag is set to True, then the the program will implement a fastq and a fasta output. IF set to false only a fasta output '''
    
    
    #Modify the outputname for the fasta file
    outputName = outputName.split('.')[0] + 'fasta'  #This converts the file extension to fasta rather than fastq. 
    
    cmd = 'usearch -fastq_filter ' + fastq + ' -fastq_maxee ' + str(maxee) + ' -fastaout ' + outputName + ' -fastq_trunclen ' + str(L)   #THe base command
    
    #if the a fastq is desired just add the addiional command to the sting 
    if fastqout == True:
        cmd = cmd + ' -fastqout ' + fastqOutputName    #The command if a fastq is desired. 
        
    #use the subprocess program to use the shell to call the uparse function.     
    sp.call(cmd.split())


def qualityFilterAndTruncateAllFastq(path, outputpath, maxee=0.5, L=400, fastqout = False):
    '''Calls usearch to quality filter reads using the maxee option and truncates the reads to a fixed length.
    If the fastq option is set to false this will also output fastq files.
    Inputs:
        Path: path to the Fastq files
        outputpath: path for the output fasta files.
        maxee: maximum expected errors. It is set to a default of 0.5
        The fastqout tag controls whether fastq files are also output. '''
    
    fastqs = [f for f in os.listdir(path) if f.endswith('fastq')]    #get the fastq file names. 

    
    #Make the directories as needed    
    try:
        os.mkdir(outputpath)
    except:
        OSError
    
    #if the fastqoutput option is set to true make a directory if needed.     
    if fastqout == True:
        fastqOutputPath = outputpath + 'Fastq'
        try:
            os.mkdir(fastqOutputPath)
        except:
            OSError
    #ecute the qualityFilterAndTrucnateAllFastqCommand.     
    for fastq in fastqs:
        qualityFilterAndTruncateFastq(join(path,fastq), join(outputpath, fastq), maxee, L, fastqout, join(fastqOutputPath, fastq))   #Calls the qulityFilterReads
    
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
 
            
def mergeAndQualityFilterFastQFiles(path, gz = False, forwardID = 'R1', reverseID ='R2', outputPath='pairedFiltered',tempPath = 'tempFasta',  truncqual = 2, minmergelen = 50, maxee = 0.5, minovlen = 10):
    """Indentifies all FastQ files present in a directory, unzipps and merges them"""
    #Inputs:
    #    path: The locations of the files
    #    outputPath: The location of the merged_Filtered Files
    #    tempFasta: the location of the unzipped and merged files, which can and should be deleted.  
    #    forwardID and reverseID: Identifier tags for the forward and reverse sequences names 
    #    gz: Tells the file if it should expect zipped gq files
    
    #Outputs: This of course makes several file transformations. It also returns the total sequence counts So we can follow how this effects my sequence counts. 
   
    #Identify all the foward and reverse files
    os.chdir(path)
    forwardFiles = [file for file in os.listdir(path) if file.find(forwardID) != -1 and file.find('fastq') != -1 and file.find('merged')==-1]   #To be placed in the first, the file must have the forward ID and also have fastq in its name. 
    reverseFiles = [file for file in os.listdir(path) if file.find(reverseID) != -1 and file.find('fastq') != -1 and file.find('merged')==-1]
    
    #Below are list for which we will add the counts of all the sequences. 
    totalCount =[]
    mergeCount = []
    filterCount = []

    #Generate the absoulte paths
    outputPath = join(path,outputPath)
    tempPath = join(path,tempPath)
#    
    #Make an output directory
    try:
        os.mkdir(join(path,outputPath))
    except:
        OSError
    
    #make temp directory for fasta files
    try: 
        
        os.mkdir(join(path, tempPath))
    except:
        OSError
    
    #Second use for loop to: decompress if needed, then merge and quality filter the files. (#With this I am assuming that the forward and reverse are in the same order). 

    
    if gz == True:
        for i in range(len(forwardFiles)):
           
            #decompress the files:
            decompressGZ(forwardFiles[i],tempPath)
            decompressGZ(reverseFiles[i],tempPath)
                        
            #Determine the proper names for all the files to be handled. 
            forward = join(tempPath,forwardFiles[i].replace('.gz', '') )
            reverse = join(tempPath, reverseFiles[i].replace('.gz','') )
            
            prefix = forwardFiles[i].split('.')[0]
            
            outputMerge = join(tempPath,prefix + '_merged.fastq' )
            outputFilter = join(outputPath, prefix + '_MF.fasta' )
            
            #Merge and filter
            #mergePairedEnds( forward, reverse, outputMerge )    The old command
            mergePairedEnds(forward, reverse, outputMerge, truncqual,  minmergelen, minovlen)         
            
            qualityFilterReads( outputMerge, outputFilter, maxee)
            #return forward, reverse, outputMerge, outputFilter
            
            #append the counts
            totalCount.append(countSeq(forward))
            mergeCount.append(countSeq(outputMerge))
            filterCount.append(countSeq(outputFilter))
       
    else:
        for i in range(len(forwardFiles)):
#            
#            #Determine the proper names for all the files to be handled. 
            forward = join(path,forwardFiles[i] )
            reverse = join(path, reverseFiles[i] )
            #print forward, reverse
            prefix = forwardFiles[i].split('.')[0]
            
            outputMerge = join(tempPath,prefix + '_merged.fastq' )
            outputFilter = join(outputPath, prefix + '_MF.fasta' )
#                    
            #Merge and filter
            #mergePairedEnds( forward, reverse, outputMerge ) old command
            mergePairedEnds(forward, reverse, outputMerge, truncqual, minmergelen, minovlen) 
            qualityFilterReads( outputMerge, outputFilter,maxee)
            
            #append the counts
            totalCount.append(countSeq(forward))
            mergeCount.append(countSeq(outputMerge))
            filterCount.append(countSeq(outputFilter))
            
    #Return the Counts for all 
    return (sum(totalCount), sum(mergeCount), sum(filterCount))
             

def countSeq(file):
    #Quick function for counting number of sequences. 
    fin = open(file)
    size = len([line for line in fin])
    fin.close
    
    if file.endswith('fastq'):
        return size/4
    if file.endswith('fasta') or file.endswith('fa'):
        return size/2  

                                                                                                                                                          
def renameReads(fasta, sampleName='Sample', index = 1):
    #This function renames the reads. It renames the reads as >Sample_index eg. >Sample_1, >Sample_2
    #Inputs:
    #       fasta = the name of the fasta file
    #       sampleName = the name of the sample prefix
    #       index =   the starting index.
    lines = [line for line in open(fasta)]
    
    
    #This goes line by line and renames the header for every sequences as desired. 
    for i in range(len(lines)):
        if lines[i].startswith('>'):
            lines[i] = '>' + sampleName + '_' + str(index) + '\n'
            index +=1
    f_out = open(fasta, 'w')
    f_out.writelines(lines)
    f_out.close()
    
    
def combineFastas(path, outputName = 'combined.fasta'):
    """Takes Fastas in an ordered List and combines them sequentially"""
    #Inputs:    
    #    path: The directory with the fasta files
    #    outputName : The name of the new combined Fasta
    
    os.chdir(path)   #need to change the path to the directory of interest.     
    files = [file for file in os.listdir(path) if file.endswith('fasta')]    #add all the fasta files to a list (shoudl be orded sequentially).
    combined = open(outputName, 'w')
        
    for file in files:    #sequentially writes the lines
        combined.writelines(open(file, 'r')) 
        combined.write('\n')      #need to tell it to skip a line. 
    combined.close()   
          
          
def renameAndCombineFastas(path, names, index = 1, outputName='combined.fasta'):
    """Calls renameReads and combineFastas to rename the sequences in each file based on the names function, and then combines them.""" 
    #Inputs:
    #   path: The path to the files
    #   names: The sample name to be applied to each sequence
    #   index: The starting number for the indexes to be passed to renameReads
    #   outputName: Obvious
    
    os.chdir(path)    #Change the dir
    files = [file for file in os.listdir(path) if file.endswith('fasta')] 
    
    #error check
    if len(names) != len(files):
        print 'The number of label names does not equal the number of Fastas'
    else:
        #rename the files
        for i in range(len(files)):
            renameReads(files[i], names[i], index)
    
        #combine the files: 
        combineFastas(path, outputName)

def dereplicateReads(fasta, output ='dereplicated.fasta'):
    #Wrapper function for Usearch that dereplicates reads. 
    #Inputs:
    #   fasta: the target merged, filtered, renamed and combined fasta
    #   output: a file name 
    
    cmd = 'usearch -derep_fulllength ' + fasta + ' -output ' + output + ' -sizeout'
    #print cmd.split()
    sp.call(cmd.split())
        
def abundanceSortReads(fasta, output = 'abundanceSorted.fasta', minSize = 2):
    #Wrapper function for Usearch that abundance sorts Reads and removes clusters that are less than a certain size. 
    #Inputs: 
    #   fasta: the target fasta to be dereplicated
    #   output: the output file name
    #   minSize: The minimum size of the clusters
 
    cmd = 'usearch -sortbysize ' + fasta + ' -output ' + output + ' -minsize ' + str(minSize)
    sp.call(cmd.split()) 

    
def clusterOtus(fasta, output ='otus.fasta', otuid =0.97):
    #Wrapper function for Usearch clusters otus and and removes chimeras.
    #Note to self, this function bassicaly removes chimeras, and then picks a rep for a cluster. It does not have any real size info that is important. 
    #Inputs: 
    #   fasta: the target fasta to be dereplicated
    #   output: the output file name
    #   The otuid sets the identity between the cluster and the sequences. Default is 0.97, it is recomended that 0.97 is minimum.
    
    cmd = 'usearch -cluster_otus ' +  fasta + ' -id ' + str(otuid)  + ' -otus ' + output
    sp.call(cmd.split())
    print cmd

    
def removeChimeras(fasta, output, database):
    #Wrapper function for Usearch that does one additional step of chimera removal
    #This function uses usearch align the otu sequencs against the reference database and determine if they are chimeras. 
    
    #Inputs: 
    #   fasta: the fasta file with the reference otu sequences. 
    #   output: the output file name
    #   database: path to the reference database for which the sequences are compared. 
    
    cmd = 'usearch -uchime_ref ' + fasta + ' -db ' + database + ' -nonchimeras '  + output + ' -strand plus' 
    sp.call(cmd.split())
    #print cmd

def mapSeqToOTUS(fasta, otuFasta, output= 'otus.uc', identity= 0.97 ):
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

def tabulateReads(path, Description= 'Standard', fout= 'numbers.txt'):
    '''This tabulates the number of reads after each run and saves the output as a tab delimted file'''
    #Note this assumes that gzip has already been handeled. 
    
    lines = ['\t'.join([Description, 'Operation','Name','Number'])] # Create a list to hold all the data The first line is the header
   
    
    #Get the names of the files. 
    forward = [f for f in os.listdir(path) if f.endswith('R1_001.fastq')]
    names = [f.split('_')[0] for f in forward]
    
    ### tabulate the reads 
    initial = [countSeqFast(join(path, f)) for f in forward]
    merged = [countSeqFast(join(path, 'tempFasta',f)) for f in os.listdir(join(path, 'tempFasta')) ]
    filtered = [countSeqFast(join(path,'pairedFiltered',f)) for f in os.listdir(join(path, 'pairedFiltered')) ]
    reads = countSeqFast(join(path,'reads.fasta'))
    otus = countSeqFast(join(path,'otus.fasta'))
    noChimera = countSeqFast(join(path,'otuFinal.fasta'))
     
    #Format the reads for writting out. Place join to the lines. We are concatenating the lists. 
    lines = lines + ['\t'.join([Description , 'Initial', names[i], str(initial[i])]) for i in range(len(initial))]    
    lines = lines + ['\t'.join([Description , 'Total Initial', 'All', str(sum(initial))])]
    
    #Merged
    lines = lines + ['\t'.join([Description , 'Merged', names[i], str(merged[i])]) for i in range(len(merged))]
    lines = lines + ['\t'.join([Description , 'Total Merged', 'All', str(sum(merged))])]
    
    #Filtered
    lines = lines + ['\t'.join([Description , 'Filtered', names[i], str(filtered[i])]) for i in range(len(filtered))]    
    lines = lines + ['\t'.join([Description , 'Total Filtered', 'All', str(sum(filtered))])]
    
    #Reads
    lines = lines + ['\t'.join([Description , 'Total Merged-Filtered', 'All', str(reads)])]
    #OTUS
    lines = lines + ['\t'.join([Description , 'OTUS', 'All', str(otus)])]
    #OTUS with out chimeras
    lines = lines + ['\t'.join([Description , 'Final OTUS', 'All', str(noChimera)])]
    
    fout = open(join(path,fout), 'w')
    fout.write('\n'.join(lines))
    fout.close()
    return lines                              

                                                                                                                  
########################################################################### Variants for doing some analysis
def fillterFastQFiles(path, outputPath='pairedFiltered',tempPath = 'tempFasta',  maxee = 0.5):
    """Indentifies all FastQ files in a directory and filters them"""
    #THis is to save time for some stuff. 
    #Inputs:
    #    path: The locations of the files
    #    outputPath: The location of the merged_Filtered Files
    #    tempFasta: the location of the unzipped and merged files, which can and should be deleted.  
    
    #    gz: Tells the file if it should expect zipped gq files
    
    #Outputs: This of course makes several file transformations. It also returns the total sequence counts So we can follow how this effects my sequence counts. 
   
    #Identify all the foward and reverse files
    os.chdir(path)
    forwardFiles = [f for f in os.listdir(path) if f.endswith('R1_001.fastq')]   
    
    #Below are list for which we will add the counts of all the sequences. 

    #Generate the absoulte paths
    outputPath = join(path,outputPath)
    tempPath = join(path,tempPath)
#    
    #Make an output directory
    try:
        os.mkdir(join(path,outputPath))
    except:
        OSError
    
    #make temp directory for fasta files
    try: 
        
        os.mkdir(join(path, tempPath))
    except:
        OSError
    
    #Second use for loop to: decompress if needed, then merge and quality filter the files. (#With this I am assuming that the forward and reverse are in the same order). 
  
    for i in range(len(forwardFiles)): 

        #print forward, reverse
        prefix = forwardFiles[i].split('.')[0]
        
        #generate the names for the files    
        Merge = join(tempPath,prefix + '_merged.fastq' )
        Filter = join(outputPath, prefix + '_MF.fasta' )
#                    
        #Merge and filter
        qualityFilterReads( Merge, Filter, maxee)
            
def mergeReads(path, tempPath = 'tempFasta',  truncqual = 5, minmergelen = 50, minovlen = 20, maxmergelen = 1000):
    """ Merges all the reads"""
    #Identify all the foward and reverse files
    os.chdir(path)
    forwardFiles = [f for f in os.listdir(path) if f.endswith('R1_001.fastq')]
    reverseFiles = [f for f in os.listdir(path) if f.endswith('R2_001.fastq')]
       
    #Generate the absoulte paths
    tempPath = join(path,tempPath)
    
    #make temp directory for fasta files
    try: 
        
        os.mkdir(join(path, tempPath))
    except:
        OSError
    
    for i in range(len(forwardFiles)):
        #Determine the proper names for all the files to be handled. 
        forward = join(path,forwardFiles[i] )
        reverse = join(path, reverseFiles[i] )
        prefix = forwardFiles[i].split('.')[0]
        
        #Make the name
        Merge = join(tempPath,prefix + '_merged.fastq' )
        #merge the reads
        
        mergePairedEnds(forward, reverse, Merge,truncqual, minmergelen )


def filterForwardReads(path, outputPath='Filtered',  maxee = 0.5):
    """Indentifies all forward FastQ files in a directory and filters them with no merging"""
    #Inputs:
    #    path: The locations of the files
    #    outputPath: The location of the merged_Filtered Files
    #    tempFasta: the location of the unzipped and merged files, which can and should be deleted.  
    #This assumes that the reads are already unzipped. 
    
    os.chdir(path)
    forwardFiles = [f for f in os.listdir(path) if f.endswith('R1_001.fastq')]   #identify the reads for filtering. 
    
    #Generate the absoulte paths
    outputPath = join(path,outputPath)
   
    #Make an output directory if required
    try:
        os.mkdir(join(path,outputPath))
    except:
        OSError

    #Second use for loop to: decompress if needed, then merge and quality filter the files. (#With this I am assuming that the forward and reverse are in the same order). 
    for i in range(len(forwardFiles)): 
        #print forward, reverse
        prefix = forwardFiles[i].split('.')[0] 
        #generate the names for the files    
        Filter = join(outputPath, prefix + '_MF.fasta' )     #this is the ouput file name   
        #FilterThe read
        qualityFilterReads( f, Filter, maxee)








                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
if __name__ == '__main__':
    
    path= '/Users/christopherrichardrivera/Documents/Pond_Project/Digester_Project/DigesterExperiment/Digester4Trimmed'
    tabulateReads(path)
    
    
    