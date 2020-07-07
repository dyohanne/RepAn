
# Including data reading and analysis funtions
source("/scratch/project_2000416/dawit/TCRanalysis/RepAnalysis26072018/RepAn/R/RepAnFns.R")
source("/scratch/project_2000416/dawit/TCRanalysis/RepAnalysis26072018/RepAn/R/RepDaAnalysisFns.R")

# Prepare CD PBMC data
CD005d0W=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/PBMC/CD005d0W.tsv")
CD005d6W=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/PBMC/CD005d6W.tsv")
CD006d0W=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/PBMC/CD006d0W.tsv")
CD006d6W=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/PBMC/CD006d6W.tsv")

CD011d0W=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/PBMC/CD011d0W.tsv")
CD011d6W=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/PBMC/CD011d6W.tsv")
CD039d0W=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/PBMC/CD039d0W.tsv")
CD039d6W=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/PBMC/CD039d6W.tsv")


# Prepare data object for analysis, and use setUp to get a repseq object
samNames=c("CD005d0W","CD006d0W","CD011d0W","CD039d0W","CD005d6W","CD006d6W","CD011d6W","CD039d6W")
CDRepseqObj <- setUp(sampleNames=samNames,samGroup=c(0,0,0,0,1,1,1,1))

devtools::use_data(CDRepseqObj, overwrite = TRUE)


# Prepare reference CDR3 dataset
referenceRepseqData <- readRDS("/scratch/project_2000416/dawit/TCRanalysis/data/referencePBMCdata.rds")

devtools::use_data(referenceRepseqData, overwrite = TRUE)


# Prepare CD Gut data

cdg1=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/Biopsy/CD1GBgfd.tsv")
cdg1T=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/Biopsy/CD1GBact.tsv")

cdg2=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/Biopsy/CD2GBgfd.tsv")
  tooShortAA <- sapply(cdg2$AMINOACID,nchar) # this sample has two AA clones which are only 2 AA long, shorter than the kmer length 3 we are working with, we just remove them, they are most likely errors.
  cdg2 <- cdg2[-which(tooShortAA < 3),]

cdg2T=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/Biopsy/CD2GBact.tsv")

cdg3=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/Biopsy/CD3GBgfd.tsv")
cdg3T=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/Biopsy/CD3GBact.tsv")

cdg4=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/Biopsy/CD4GBgfd.tsv")
cdg4T=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/Biopsy/CD4GBact.tsv")

cdg5=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/Biopsy/CD5GBgfd.tsv")
cdg5T=readSample("/scratch/project_2000416/dawit/TCRanalysis/data/Biopsy/CD5GBact.tsv")

# Prepare data object for analysis, and use setUp to get a repseq object

samNames=c("cdg1","cdg2","cdg3","cdg4","cdg5","cdg1T","cdg2T","cdg3T","cdg4T","cdg5T")

CDGutRepseqObj <- setUp(sampleNames=samNames,samGroup=c(0,0,0,0,0,1,1,1,1,1))

devtools::use_data(CDGutRepseqObj, overwrite = TRUE)

