\name{extractDASubRepertoire}
\alias{extractDASubRepertoire}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractDASubRepertoire(repSeqObj, cutoff = 0.1, method = "pvalue")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repSeqObj}{
%%     ~~Describe \code{repSeqObj} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (repSeqObj, cutoff = 0.1, method = "pvalue") 
{
    if (method == "pvalue") {
        colToCheck = "pval"
    }
    else {
        colToCheck = "padj"
    }
    significantSubRepertoires = rownames(repSeqObj$cDaResult)[which(repSeqObj$cDaResult[, 
        colToCheck] < cutoff)]
    if (!length(significantSubRepertoires) >= 1) {
        cat("No differentially abundant subrepertoires found with given cutoff!\n")
        return(repSeqObj)
    }
    DAClonotypes = NULL
    for (subR in significantSubRepertoires) {
        sigSubRepIndex <- repSeqObj$clusterMatchTable[rownames(repSeqObj$clusterMatchTable) == 
            subR, ]
        sigSubRepIndex <- sigSubRepIndex[!is.na(sigSubRepIndex)]
        for (s in names(sigSubRepIndex)) {
            sampleSelected = repSeqObj$scaledSampleData[[s]][repSeqObj$withinSampleClusters[[s]]$clusters == 
                sigSubRepIndex[s], ]
            datemp = cbind(sampleSelected$AMINOACID, rep(subR, 
                length(sampleSelected$AMINOACID)))
            DAClonotypes = rbind(DAClonotypes, datemp)
        }
    }
    DAClonotypeAbundanceMatrix = NULL
    for (s in names(repSeqObj$scaledSampleData)) {
        DAClonotypeAbundanceMatrix = cbind(DAClonotypeAbundanceMatrix, 
            repSeqObj$scaledSampleData[[s]][match(DAClonotypes[, 
                1], repSeqObj$scaledSampleData[[s]]$AMINOACID), 
                c("COUNT")])
    }
    colnames(DAClonotypes) <- c("AMINOACID", "Subrepertoire")
    rownames(DAClonotypeAbundanceMatrix) = DAClonotypes[, 1]
    colnames(DAClonotypeAbundanceMatrix) = names(repSeqObj$scaledSampleData)
    DAClonotypeAbundanceMatrix[is.na(DAClonotypeAbundanceMatrix)] <- 0
    repSeqObj = addItemToObject(repSeqObj, significantSubRepertoires, 
        "cDaClusters")
    repSeqObj = addItemToObject(repSeqObj, DAClonotypes, "cDaClonotypesWithCluster")
    repSeqObj = addItemToObject(repSeqObj, DAClonotypes[, 1], 
        "cDaClonotypes")
    repSeqObj = addItemToObject(repSeqObj, DAClonotypeAbundanceMatrix, 
        "cDaClonotypeAbundanceMatrix")
    return(repSeqObj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
