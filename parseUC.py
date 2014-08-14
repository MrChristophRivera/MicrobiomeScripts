###Scripts for parsing the Uparse information to be compatible with the Qiime
import os


#note to self that we need to strip the white space, and then split based on tab

def generateOtuDict(ucFile):
    """Generates a dictionary containing the otu as key and the sequence names as list"""
    lines = [line.strip() for line in open(ucFile)]    #place all the lines into a list and remove the white space
    
    #Create the dictionary
    otuD = {}
    
    matched = 0.0   #number to determine how many matched

    #parse out the pertinant information and place into the dictionary
    for line in lines:
        line = line.split('\t')
        if line[9]!= '*':
            otuD.setdefault(line[9],[]).append(line[8])
            matched +=1

    print    'Fraction Matched: ', matched/len(lines)   #print out the number matched.
    return otuD

def makeOtuList(otuD):
    """Takes the an OtuDictionary and colates it into a tuppled list that can be sorted. """
    
    #create the tuppled list and sort it. 
    otuList = [ (int(key.split('_')[1]), key, otuD[key]) for key in otuD]   #This is sorted based on the key value int. 
    otuList.sort()
    return otuList
    

def generateMapLine(otu):
    """ Configures the data from the Otu List into a writable format"""
    """Creates a list with the the denovo entry first, and then the sequence name."""
    #The list is then combined into a string with tab seperation. 
    MapLine = [otu[1]] + otu[2]
    MapLine = '\t'.join(MapLine) +'\n'
    return MapLine
    
def writeOtus(otuList, outputName= 'seqs_otus.txt'):
    """Write the Otus in the MapLines to a tab delimited text file.""" 
    fin = open(outputName, 'w')
    for otu in otuList:
        line = generateMapLine(otu)
        fin.write(line)
    
    fin.close()
 
def parseUCandWriteOtus(ucFile, outputName='seqs_otus.txt'):
    """This client function parses a UC file and writes otu a tab delimted otufile"""
    
    otuD = generateOtuDict(ucFile)  
    otuList = makeOtuList(otuD)
    writeOtus(otuList, outputName)
    
    
    
    
    
    
        