# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 15:11:23 2014

@author: chrisrichardrivera
"""


from os.path import join, split
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from os.path import join
from os import getcwd, chdir, listdir

def convertQscores(qual):
    '''Converts the quality score to an number.
    Parameters: 
        qual (str): The quality scores
    Returns:
        list with numeric qscores'''
    return [ord(q)-33 for q in qual]

def getNumericQscoresFromFastq(fastq):
    '''Gets all the qscores from a fast q file and converts them to numeric.
    Parameters:
        fastq (str): Path to fastq file.
    Returns:
        list of list'''
    with open(fastq) as f:
        return [convertQscores(qual) for (title, seq, qual) in  FastqGeneralIterator(f)]
    


if __name__=='__main__':
    
    path = '/Users/chrisrichardrivera/Desktop/Digester1Trimmed/Sequences'
    fastqs = [join(path, f) for f in listdir(path) if f.endswith('fastq')]

    qs = getNumericQscoresFromFastq(fastqs[0])
            