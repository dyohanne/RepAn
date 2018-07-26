# RepAn
Analysis functions for deep sequenced immune repertoire (Repseq) datasets in R

RepAn mainly performs differential abundance analysis to identify condition associated TCR CDR3beta sequences by comparing samples 
in two treatment/condition groups. It uses unsupervised within sample CDR3 clustering to first define clusters of highly similar
CDR3 sequences, then finds their closest matches across samples. Once matching clusters of CDR3s across samples are found, which are
referred to as sub-repertoires, it applies statistcal testing to identify diffferentially abundant (DA) sub-repertoires by comparing
the abundances of sub-repertoires between condition groups. CDR3 sequences belonging to DA sub-repertoires are collected and 
put through ranking based filtering to identify most likely condition associated CDR3 sequences. 

For detailed description of the differential abundance analysis method, please see (link)

To install RepAn package in R, first install R devtools package and run the following in R: devtools::install_github("dyohanne/RepAn")


RepAn for now accepts only [immunoseq](https://www.adaptivebiotech.com/immunoseq) formated datasets.

To perform differential abundance analysis, first read in the repertoire samples using the readSample function as follows:

```
cd1=readSample("data/PBMC/CD005d0W.tsv")
cd1T=readSample("data/PBMC/CD005d6W.tsv")
cd2=readSample("data/PBMC/CD006d0W.tsv")
cd2T=readSample("data/PBMC/CD006d6W.tsv")

cd3=readSample("data/PBMC/CD011d0W.tsv")
cd3T=readSample("data/PBMC/CD011d6W.tsv")
cd4=readSample("data/PBMC/CD039d0W.tsv")
cd4T=readSample("data/PBMC/CD039d6W.tsv")
```

```
CXX=clang++
```



RepAn also implements basic functions for descriptive function for the analysis of immune repertoire datasets. These enable:
- basic repertoire summary and statistics analysis 
- repertoire overlap analysis
- repertoire diversity analysis


