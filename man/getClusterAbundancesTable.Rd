\name{getClusterAbundancesTable}
\alias{getClusterAbundancesTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getClusterAbundancesTable(repSeqObj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repSeqObj}{
%%     ~~Describe \code{repSeqObj} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (repSeqObj) 
{
    if ("clusterMatchTable" \%in\% names(repSeqObj)) {
        clusMatchTable = repSeqObj$clusterMatchTable
    }
    else {
        stop("Cluster match table not found !.")
    }
    cAbundanceTable = NULL
    cNumberOfClones = NULL
    cRelativeCloneSizeFrequencyTable = NULL
    cRelativeAverageCloneSizeTable = NULL
    for (subRep in 1:nrow(clusMatchTable)) {
        subRepIndex <- clusMatchTable[subRep, ]
        subRepertoireAbundance = c()
        subRepertoireNumberOfClonotypes = c()
        subRepertoireRelativeCloneSizeFreq = c()
        subRepertoireRelativeAverageCloneSize = c()
        for (s in names(subRepIndex)) {
            sampleSelected = repSeqObj$scaledSampleData[[s]][repSeqObj$withinSampleClusters[[s]]$clusters == 
                subRepIndex[s], ]
            subRepertoireAbundance = c(subRepertoireAbundance, 
                sum(sampleSelected$COUNT))
            if (is.na(subRepIndex[s])) {
                subRepertoireNumberOfClonotypes = c(subRepertoireNumberOfClonotypes, 
                  NA)
            }
            else {
                subRepertoireNumberOfClonotypes = c(subRepertoireNumberOfClonotypes, 
                  nrow(sampleSelected))
            }
            subRepertoireRelativeCloneSizeFreq = c(subRepertoireRelativeCloneSizeFreq, 
                sum(sampleSelected$COUNT)/sum(repSeqObj$scaledSampleData[[s]]$COUNT))
            samSelectedCountPerClone = sum(sampleSelected$COUNT)/nrow(sampleSelected)
            samCountPerClone = sum(repSeqObj$scaledSampleData[[s]]$COUNT)/nrow(repSeqObj$scaledSampleData[[s]])
            samSelectedVsTotal = samSelectedCountPerClone/samCountPerClone
            subRepertoireRelativeAverageCloneSize = c(subRepertoireRelativeAverageCloneSize, 
                samSelectedVsTotal)
        }
        cAbundanceTable = rbind(cAbundanceTable, subRepertoireAbundance)
        cNumberOfClones = rbind(cNumberOfClones, subRepertoireNumberOfClonotypes)
        cRelativeCloneSizeFrequencyTable = rbind(cRelativeCloneSizeFrequencyTable, 
            subRepertoireRelativeCloneSizeFreq)
        cRelativeAverageCloneSizeTable = rbind(cRelativeAverageCloneSizeTable, 
            subRepertoireRelativeAverageCloneSize)
    }
    colnames(cAbundanceTable) = colnames(clusMatchTable)
    rownames(cAbundanceTable) = 1:nrow(clusMatchTable)
    repSeqObj <- addItemToObject(repSeqObj, cAbundanceTable, 
        "cAbundanceTable")
    colnames(cNumberOfClones) = colnames(clusMatchTable)
    rownames(cNumberOfClones) = 1:nrow(clusMatchTable)
    repSeqObj <- addItemToObject(repSeqObj, cNumberOfClones, 
        "cNumberOfClones")
    colnames(cRelativeCloneSizeFrequencyTable) = colnames(clusMatchTable)
    rownames(cRelativeCloneSizeFrequencyTable) = 1:nrow(clusMatchTable)
    repSeqObj <- addItemToObject(repSeqObj, cRelativeCloneSizeFrequencyTable, 
        "cRelativeAbundanceTable")
    colnames(cRelativeAverageCloneSizeTable) = colnames(clusMatchTable)
    rownames(cRelativeAverageCloneSizeTable) = 1:nrow(clusMatchTable)
    repSeqObj <- addItemToObject(repSeqObj, cRelativeAverageCloneSizeTable, 
        "cRelativeAverageCloneSizeTable")
    return(repSeqObj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
