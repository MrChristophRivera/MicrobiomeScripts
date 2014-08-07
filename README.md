MicrobiomeScripts
=================
These are scripts for interfacing the Uparse algorithms with Python and Qiime for 16s Microbiome Analysis. 

There are 5 files:

UparseWrappers contains wrapper functions that call the Uparse commands via a python interface. 

UparsePipeline calls the Wrapper functions for calling the Uparse Wrapper functions. 

parseUC contains functions for parse the UC file into a tab delimited file that can be fed into the qiime pathway to make an otu table. 

Qiime Wrappers contains wrapper functions for calling Qiime easily via ipython notebook. 

filterFastq contains various functions for filtering and plotting fastq files.
