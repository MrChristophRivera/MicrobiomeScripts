#Scripts to parse filter reads. 

from Bio import SeqIO
import os
from os.path import join
import numpy as np
import matplotlib.pyplot as mp
#from numpy.random import randint
import gzip

def getQscore(seq):
    #This script returns the Qscore
    return seq.letter_annotations["phred_quality"]

def testQscores(subQ, q =3):
    #This function determines if the all of q scores within in a subsequnce are greater than the threshold. Returns True if True otherwise returns false. 
    #Note that it does not test the first bp (which is good since it is typically low quality (we can consider filtering this one out). 
    
    subQ = np.array(subQ)  #convert to an array
    
    return sum(subQ >= q) == len(subQ) #If every subq score is greater than the threshold return True otherwise return false. 

def determineLastHighQualityBP(seq, r = 3,q = 3):
    #This function determines the position of the last usable base pair based on a moving window of q scores. 
    #It moves along the qscore and determines if the next r sequences are higher than the threshold, if they are higher than the q score threshold it advances. 
    #If not it returns the value of that index
    
    #use a while loop to go through the sequence and test the function
    #use testQscores to test subsequence, if the next subsequence passes, advance otherwise return the index
    
    Qscores=getQscore(seq)  #call getQscores to get the Qscores
    
    i = 0
    while i < len(Qscores)-r:
        subQ = Qscores[i+1:i+r+1]   #Index the subQ scores spanning from i+1 to i+r
        if testQscores(subQ, q):   #if it the Q scores pass advance, if not return the index + 1 (this is to make sure that the sequence includes proper sequence). 
            i += 1
        else:
            return i +1     
    return len(seq)  #if the function reaches the end of the sequence up to r, return the final index. 
   
def testAlgorithm(Qscores, r=3, q=3):
    #This is allowed for me to 
    i = 0
    while i < len(Qscores)-r:
        subQ = Qscores[i+1:i+r+1]   #Index the subQ scores spanning from i+1 to i+r
        if testQscores(subQ, q):   #if it the Q scores pass advance, if not return the index + 1 (this is to make sure that the sequence includes proper sequence). 
            i += 1
        else:
            return i +1       
    return len(Qscores) 
    

def trimSeq(seq, r =3, q = 3):
    #This script trims a sequence based on its fastq score. 
    #It follows the qscore via a moving window #starting at bp 1. 
    #If the next r errors probabilities (Q score) are lower than the qualtity threshold, 
    #it truncates the read at that base. 
    #Inputs are the sequences (fastq format), 
    #r (which is the run length parameter) (default of 3) 
    #and the Q score threshold (q) with default of 3 
    
    lastBP = determineLastHighQualityBP(seq,r,q)     #determine the last HighQuality bp
    return seq[0:lastBP]    #Return the properly truncated sequence. 

def greaterThanMinLength(seqLength, trimmedSeqLength, p=0.75):
    #simple function determines if the trimed sequence lenth is greater than the min sequence length. 
    #seqLength is the sequence length, index is the index returned by the determineLastHighQualityBP, p is minimum percentage. 
    return (trimmedSeqLength >= p*seqLength)
    
def lessThanMaxN(trimmedSeq, maxN = 0):
    #Returns True if the sequences has less or or equal to max N.  
    s = trimmedSeq.seq    #get the sequence as a string
    return (s.count('N') <= maxN)
                
def filterSeq(seq,trimmedSeq, p = 0.75, maxN = 0):
    #this uses the above scripts to determine whether the sequence should be fitlered or not. 
    
    #determine the sequence lenths
    seqLength= len(seq)
    trimmedSeqLength = len(trimmedSeq)
    return (greaterThanMinLength(seqLength, trimmedSeqLength, p) & lessThanMaxN(trimmedSeq, maxN))
  
def trimAndFilterSeqs(seqs, r=3, q=3, p=0.75, maxN = 0):
    #This script trims and fitlers sequences using the above functions and appends them to a new list. 
    #if the sequences is filtered, then the will be appended as a shorted sequences, that should be pruned out during pandaseq
    seqsTrimmedAndFiltered=[]
    count = 0
    for seq in seqs:
        seqT = trimSeq(seq,r,q)
        if filterSeq(seq, seqT, p,maxN):
            seqsTrimmedAndFiltered.append(seqT)
            count+=1
        else:
            seqsTrimmedAndFiltered.append(seq[0:1])     #if this does not pass the filtering process, just append the first one. #We will require an overlap of >1 in the pandaseq step/ 
            
    return seqsTrimmedAndFiltered, count
        


def trimAndFilterFiles(directory, outputDir = 'Filtered', compress = True,  r=3, q=3, p=0.75, maxN=0):
    #This script trims and fitlers all the sequences with in a list of files using the above scripts. 
    #it assumes that the files are gzipped fastq files. ( i can modify it as needed). 
    #if compress is True, then use gzip to compress
    
    #Determine the files, access them using gzip then parse them using Biopython and the above functions. 
    
    files = [file for file in os.listdir(directory) if file.endswith('fastq.gz')]
    try:
        os.mkdir(outputDir)
    except OSError:
        print'The directory exists!'
        
    
    sequencesPassed = 0.0
    totalSequences = 0.0
    
    #go file by file and trim. 
    #Write the sequences- if the sequences fail the trim then they will be truncated. 
    
    for file in files:
        #print file            
        f =  gzip.open(join(directory,file))     #create a handle for the file and access with gz
        
        seqs = [seq for seq in SeqIO.parse(f, 'fastq')]
        totalSequences += len(seqs)
        seqs , count = trimAndFilterSeqs(seqs, r, q, p, maxN)
        sequencesPassed += count
        
        #Write the sequences to a new file. Compress the files if desired above. (I likely will be doing this)
        if compress:        
            h = gzip.open(join(directory,'Filtered',file), 'w')
            SeqIO.write(seqs,h,'fastq')
            h.close()
        else: 
            SeqIO.write(seqs, join(directory, 'Filtered', file.replace('.gz', '')), 'fastq')
        f.close()    #close the file to conserve memory. 
    
    print 'Total Sequences:  ', totalSequences
    print 'SequencesPassed:  ', sequencesPassed
    print sequencesPassed/totalSequences*100
        
  
def getQscoresFromFiles(directory):
    #This function acesses files and get qscores for plotting (this is for fun really) #it assumes that the fastq files are compresed
    
    files = [file for file in os.listdir(directory) if file.endswith('fastq.gz') or file.endswith('fastq')]
    
    
    Qscores =[]
    
    #append the qscores
    
    for file in files:
        #determine the proper way to open the file: if its compressed, use gzip. if not use the normal compression. 
        
        if file.endswith('fastq.gz'):
            file = gzip.open(join(directory,file))
        else: 
            file = open(join(directory,file))
        for seq in SeqIO.parse(file,'fastq'):
            Qscores.append(getQscore(seq))
            
    
    return Qscores

def createQscoreDict(Qscores):
    #Go through a list of Qscores and determine make a dictionary containing the the qscores for every bp. 
    Qdict ={}
    for qscore in Qscores:
        for i in range(len(qscore)):
            Qdict.setdefault(i+1,[]).append(qscore[i])
    
    return Qdict
        
def plotMeanQscore(Qdict):
      #this function plots the mean values of a qscore
      meanQscores = []
      for key in Qdict:
          Mean = np.mean(Qdict[key])
          StDev = np.std(Qdict[key])
          meanQscores.append((key, Mean, StDev))
      meanQscores.sort()
      
      #plot the qscores
      bp,Mean, StDev = zip(*meanQscores)
      mp.errorbar(bp, Mean, yerr=StDev, 
         fmt='-', ecolor='black', elinewidth=0.5, capsize=3,
         barsabove=False, lolims=False, uplims=False,
         xlolims=False, xuplims=False, errorevery=1,
         capthick=None ,linewidth=2.5, color ='green')
      mp.xlabel('Base Pair', fontsize=16)
      mp.ylabel('Mean Qscore',fontsize=16)
      mp.title('Mean Qscore versus Base Pair',fontsize=20)
      mp.show()
      
      
def plotMeanQscoreLength(Qscores):
    #This plots the mean qscore length as a histogram, which is really just a proxy for sequence lenth
    qlengths = [len(qs) for qs in Qscores]
    mp.hist(qlengths, bins = np.max(qlengths)/2)
    mp.xlabel('Sequence Length',fontsize=16)
    mp.ylabel('Number', fontsize=16)
    mp.show()
      
def aggregateQscores(qscores):
    #this function just plots a histogram of the qscores.
    Qscores = []
    for q in qscores:
           Qscores +=q
    return Qscores


def histogramQscoresBulk(qscores):
    #this draws a histgram of the qscores... for curiosity sake. Of bulk.
    mp.hist(qscores, bins = 40)
    mp.xlabel('Qscore',fontsize=16)
    mp.ylabel('Number', fontsize=16)
    mp.show()

def getQscoresFromFile(file, directory):
    #Gets the Qscores from one indivial file for comparison. 
    
    Qscores =[]   #make a qscores list
    
    #Create a file handle as needed
    if file.endswith('fastq.gz'):
        file = gzip.open(join(directory,file))
    else: 
        file = open(join(directory,file))
        
    #Parse out the Qscores.     
    for seq in SeqIO.parse(file,'fastq'):
        Qscores.append(getQscore(seq))
            
    
    return Qscores
    

if __name__== "__main__":
    
    directory = '/Users/christopherrichardrivera/Documents/Pond_Project/Digester_Project/Digester-November_11-2013/SequencingData/pairedFiltered'

    
    Qscores = getQscoresFromFiles(directory)
    
    
    Qdict = createQscoreDict(Qscores)
    
    plotMeanQscore(Qdict)
    mp.figure()
    plotMeanQscoreLength(Qscores)
    

def testSeqIO(file, directory):
    qscores =[]
    seqs = [seq for seq in SeqIO.parse(join(directory, file), 'fastq')]
    for seq in seqs:
        qscores.append(getQscore(seq))
    
    print "Number of Qscores: ", len(qscores)
    print  "Number of sequences: ", len(seqs)
    
    
    
    