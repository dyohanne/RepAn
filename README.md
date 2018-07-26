# RepAn
Analysis functions for deep sequenced immune repertoire (Repseq) datasets in R

RepAn mainly performs differential abundance analysis to identify condition associated TCR CDR3beta sequences by comparing samples 
in two treatment/condition groups. It uses unsupervised within sample CDR3 clustering to first define clusters of highly similar
CDR3 sequences, then finds their closest matches across samples. Once matching clusters of CDR3s across samples are found, which are
referred to as sub-repertoires, it applies statistcal testing to identify diffferentially abundant (DA) sub-repertoires by comparing
the abundances of sub-repertoires between condition groups. CDR3 sequences belonging to DA sub-repertoires are collected and 
put through ranking based filtering to identify most likely condition associated CDR3 sequences. 

For detailed description please see (link)

RepAn also implements basic functions for descriptive function for the analysis of immune repertoire datasets. These enable:
- basic repertoire summary and statistics analysis 
- repertoire overlap analysis
- repertoire diversity analysis

RepAn for now accepts only [immunoseq](https://www.adaptivebiotech.com/immunoseq) format datasets.

