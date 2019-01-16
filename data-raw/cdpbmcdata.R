
# Including data reading and analysis funtions
source("/wrk/dawit/TCRanalysis/RepAnalysis26072018/RepAn/R/RepAnFns.R")
source("/wrk/dawit/TCRanalysis/RepAnalysis26072018/RepAn/R/RepDaAnalysisFns.R")

# read in data
CD005d0W=readSample("/wrk/dawit/TCRanalysis/data/PBMC/CD005d0W.tsv")
CD005d6W=readSample("/wrk/dawit/TCRanalysis/data/PBMC/CD005d6W.tsv")
CD006d0W=readSample("/wrk/dawit/TCRanalysis/data/PBMC/CD006d0W.tsv")
CD006d6W=readSample("/wrk/dawit/TCRanalysis/data/PBMC/CD006d6W.tsv")

CD011d0W=readSample("/wrk/dawit/TCRanalysis/data/PBMC/CD011d0W.tsv")
CD011d6W=readSample("/wrk/dawit/TCRanalysis/data/PBMC/CD011d6W.tsv")
CD039d0W=readSample("/wrk/dawit/TCRanalysis/data/PBMC/CD039d0W.tsv")
CD039d6W=readSample("/wrk/dawit/TCRanalysis/data/PBMC/CD039d6W.tsv")


# Prepare data object for analysis, and use setUp to get a repseq object
samNames=c("CD005d0W","CD006d0W","CD011d0W","CD039d0W","CD005d6W","CD006d6W","CD011d6W","CD039d6W")
CDRepseqObj <- setUp(sampleNames=samNames,samGroup=c(0,0,0,0,1,1,1,1))

devtools::use_data(CDRepseqObj, overwrite = TRUE)


# Prepare reference CDR3 dataset
referenceRepseqData <- readRDS("/wrk/dawit/TCRanalysis/data/referencePBMCdata.rds")

devtools::use_data(referenceRepseqData, overwrite = TRUE)
