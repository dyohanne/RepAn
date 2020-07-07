# RepAn
RepAn is an R package that implements a differential abundance analysis method for deep sequenced TCR immune repertoire (Repseq) datasets to identify enriched/expanded clonotypes associated with a condition/disease.

RepAn mainly identifies condition associated TCR CDR3beta sequences by comparing samples 
in two treatment/condition groups. It uses unsupervised within sample CDR3 clustering to first define clusters of highly similar
CDR3 sequences, then finds their closest matches across samples. Once matching clusters of CDR3s across samples are found, which are
referred to as sub-repertoires, it applies statistcal testing to identify diffferentially abundant (DA) sub-repertoires by comparing
the abundances of sub-repertoires between condition groups. CDR3 sequences belonging to DA sub-repertoires are collected and 
put through a ranking based filtering to identify the most likely condition associated CDR3 sequences using permutation to calculate statistical significance of the observed rankings. 

For detailed description of the differential abundance analysis method, please see (https://www.biorxiv.org/content/early/2018/12/07/490102)

## Installation
To install RepAn package in R, first install R devtools package and run its install_github function as follows: 

```
install.packages("devtools")
devtools::install_github("dyohanne/RepAn")
```
If installing R package dependencies with the previous command gives errors, first install RepAn without the dependencies and use the package installation functions that come with RepAn as follows:

```
devtools::install_github("dyohanne/RepAn",dependencies=F)
loadPacks()
loadBioconductorPacks()
```
We recommend the installation of parallelRandomForest to speed up the RandomForest based clonotype ranking task in RepAn, especially for big datasets. Please follow the instruction on (https://bitbucket.org/mkuhn/parallelrandomforest/src/parallelRandomForest/) to install parallelRandomForest or download the source of parallelRandomForest and install locally in R.    

## Running RepAn
RepAn is mainly designed to work with genomic TCR repertoire datasets and accepts [immunoseq](https://www.adaptivebiotech.com/immunoseq) formatted datasets. It also accepts MiXCR format data, in which case, if the data comes from cDNA libraries, the analysis can be interpreted as differential expression of TCRs.

To perform differential abundance analysis, first read in the repertoire samples using the readSample function (for immunoseq format data) or the readMiXCR function (for MiXCR format data). For example to read condition 1 samples named sample1.tsv, sample2.tsv,... , and condition 2 samples named sample1T.tsv, ..., do the following ; 

```
s1=readSample("sample1.tsv")
s2=readSample("sample2.tsv")
s3=readSample("sample3.tsv")

s1T=readSample("sample1T.tsv")
s2T=readSample("sample2T.tsv")
s3T=readSample("sample3T.tsv")
```

Then prepare the vector of sample names, and run the Repseq data set up function by passing the sample names vector and the group or condition status vector (0s for condition 1,and 1s for condition 2 samples). The setUp function returns a Repseq data object that contains all the samples normalized to counts per million.

```
samNames=c("s1","s2","s3","s1T","s2T","s3T")

repObj <- setUp(sampleNames=samNames,samGroup=c(0,0,0,1,1,1))
```

Now we can run differential abundance analysis by passing the Repseq data object to the runDaAnalysis function along with other parameters. Importantly,clusterby="NT" to use nt based 4-mer frequency to represent CDR3 sequences, resampleSize to define the downsampling size (default 5000), repeatResample=T to perform repeat resampling in the analysis, nRepeats the number of repeat resamples, etc. Read the help for runDaAnalysis function for details. The function returns a list with three elements unless returnAll is set to false: the first two elements are data frames of all candidate DA CDR3s along with their abundance in the samples and DA filtering statistics from enrichment and de-enrichment analyses respectively, third element is the address of the intermediate results for the repeat resample runs of the anaysis. 

```
results <- runDaAnalysis(repObj,clusterby="NT",kmerWidth=4,paired=T,clusterDaPcutoff=0.1,positionWt = F,distMethod="euclidean",matchingMethod="km",nRepeats=10,resampleSize=5000,useProb=T,returnAll=T,nRR=1000)
```

Repseq data object called CDRepseqObj prepared from PBMC TCRB CDR3 repertoire datasets of four celiac disease (CD) patients before and after a 3 day oral gluten challenge is included in RepAn. RepAn also includes data CDGutRepseqObj from CD Gut dataset (n=5) thas has gut biopsy TCR repertoires of celiac disease patients during active disease and after one year long gluten-free diet. We can use CDRepseqObj to find CD associated CDR3 sequences enriched in the repertoires following a 3-day gluten challenge as follows:

```
# Check the members of the CDRepseqObj object
names(CDRepseqObj)

# Names of pre-gluten challenge day 0 samples
CDRepseqObj$samNames[CDRepseqObj$group==0]

# Names of post-gluten challenge day 6 samples
CDRepseqObj$samNames[CDRepseqObj$group==1]

# We finally run RepAn DA analysis function on the data with nt 4-mer based clustering, 5 analysis with downsampled repertoire size of 3000 per sample 
results <- runDaAnalysis(CDRepseqObj,clusterby="NT",kmerWidth=4,paired=T,clusterDaPcutoff=0.1,positionWt=F,
  distMethod="euclidean",matchingMethod="km",nRepeats=5,resampleSize=3000,useProb=T,returnAll=T,nRR=1000)
```

Significantly differentially enriched CDR3s (with permuted p-values < 0.05 and q-value < 0.05), that have shown significant clonal expansion in condition 2 (gluten challenged samples), can be accessed from the result using :

```
enrichedCDR3s <- TopDAClonotypes(results,enriched=T,pValueCutoff=0.05,qValue=0.05)

```

Significantly differentially de-enriched CDR3s (with permuted p-values < 0.05 and q-value < 0.05), that have shown significant clonal contraction in condition 2, can be accessed from the result using :

```
de_enrichedCDR3s <- TopDAClonotypes(results,enriched=F,pValueCutoff=0.05,qValue=0.05)
```



RepAn also contains basic functions for descriptive analysis of immune repertoire datasets that enable:
- basic repertoire summary and statistical analysis 
- repertoire overlap analysis
- repertoire diversity analysis
- V- and J- gene usage analysis



