# RepAn
Analysis functions for deep sequenced immune repertoire (Repseq) datasets in R

RepAn mainly performs differential abundance analysis to identify condition associated TCR CDR3beta sequences by comparing samples 
in two treatment/condition groups. It uses unsupervised within sample CDR3 clustering to first define clusters of highly similar
CDR3 sequences, then finds their closest matches across samples. Once matching clusters of CDR3s across samples are found, which are
referred to as sub-repertoires, it applies statistcal testing to identify diffferentially abundant (DA) sub-repertoires by comparing
the abundances of sub-repertoires between condition groups. CDR3 sequences belonging to DA sub-repertoires are collected and 
put through ranking based filtering to identify most likely condition associated CDR3 sequences. 

For detailed description of the differential abundance analysis method, please see (link)

To install RepAn package in R, first install R devtools package and run its install_github function as follows: 

```
install.packages("devtools")
devtools::install_github("dyohanne/RepAn")
```

RepAn for now accepts only [immunoseq](https://www.adaptivebiotech.com/immunoseq) formated datasets.

To perform differential abundance analysis, first read in the repertoire samples using the readSample function. For example to read condition 1 samples named sample1.tsv, ... , and condition 2 samples named sample1T.tsv, ..., do the following ; 

```
s1=readSample("sample1.tsv")
s2=readSample("sample2.tsv")
s3=readSample("sample3.tsv")

s1T=readSample("sample1T.tsv")
s2T=readSample("sample2T.tsv")
s3T=readSample("sample3T.tsv")
```

Then prepare the vector of sample names, and run the Repseq data set up function by passing the sample name vector and the group status vector (0s for condition 1,and 1s for condition two samples). The setUp function returns a Repseq data object that contains all the samples normalized to counts per million.

```
samNames=c("s1","s2","s3","s1T","s2T","s3T")

repObj <- setUp(sampleNames=samNames,samGroup=c(0,0,0,1,1,1))
```

Now we can run differential abundance analysis by passing the Repseq data object to the runDaAnalysis function along with other parameters. Importantly,clusterby="NT" to use nt based 4-mer frequency to represent CDR3 sequences, resampleSize to define the downsampling size (default 5000), repeatResample=T to perform repeat resampling in the analysis, nRepeats the number of repeat resamples etc. Read the help for runDaAnalysis function for details. The function returns a list with two elements unless returnAll is set to false: the first element is a data frame of all candidate DA CDR3s along with their abundance in the samples and DA filtering statistics, second element is a list of intermediate results from each repeat resample run of the anaysis. 

```
results <- runDaAnalysis(repObj,clusterby="NT",kmerWidth=4,paired=T,clusterDaPcutoff=0.1,positionWt = F,distMethod="euclidean",matchingMethod="km",nRepeats=10,resampleSize=5000,useProb=T,returnAll=T,nRR=1000)
```

Significantly differentially enriched CDR3s (with permuted p-values < 0.05), that have shown significant clonal expansion in condition 2, can be accessed from the result using :

```
enriched <- results[[1]][results[[1]]$permutedEnPval < 0.05,]
```

Significantly differentially de-enriched CDR3s (with permuted p-values < 0.05), that have shown significant clonal contraction in condition 2, can be accessed from the result using :

```
de-enriched <- results[[1]][results[[1]]permutedDeEnPval < 0.05,]
```


RepAn also implements basic functions for descriptive analysis of immune repertoire datasets that enable:
- basic repertoire summary and statistics analysis 
- repertoire overlap analysis
- repertoire diversity analysis
- etc.


