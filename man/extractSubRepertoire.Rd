\name{extractSubRepertoire}
\alias{extractSubRepertoire}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extractSubRepertoire(repSeqObj, subReps = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repSeqObj}{
%%     ~~Describe \code{repSeqObj} here~~
}
  \item{subReps}{
%%     ~~Describe \code{subReps} here~~
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
function (repSeqObj, subReps = NULL) 
{
    if (is.null(subReps)) 
        stop("No subrepertoires given.")
    if (!item.exists(repSeqObj, "clusterMatchTable")) {
        stop("clusterMatchTable was not detected.")
    }
    subRepClonotypes = c()
    for (subR in subReps) {
        if (!subR \%in\% rownames(repSeqObj$clusterMatchTable)) 
            next
        sigSubRepIndex <- repSeqObj$clusterMatchTable[rownames(repSeqObj$clusterMatchTable) == 
            subR, ]
        sigSubRepIndex <- sigSubRepIndex[!is.na(sigSubRepIndex)]
        for (s in names(sigSubRepIndex)) {
            sampleSelected = repSeqObj$scaledSampleData[[s]][repSeqObj$withinSampleClusters[[s]]$clusters == 
                sigSubRepIndex[s], ]
            subRepClonotypes = c(subRepClonotypes, sampleSelected$AMINOACID)
        }
    }
    DAClonotypeAbundanceMatrix = NULL
    for (s in names(repSeqObj$scaledSampleData)) {
        DAClonotypeAbundanceMatrix = cbind(DAClonotypeAbundanceMatrix, 
            repSeqObj$scaledSampleData[[s]][match(subRepClonotypes, 
                repSeqObj$scaledSampleData[[s]]$AMINOACID), c("COUNT")])
    }
    rownames(DAClonotypeAbundanceMatrix) = subRepClonotypes
    colnames(DAClonotypeAbundanceMatrix) = names(repSeqObj$scaledSampleData)
    DAClonotypeAbundanceMatrix[is.na(DAClonotypeAbundanceMatrix)] <- 0
    repSeqObj = addItemToObject(repSeqObj, subRepClonotypes, 
        "cSubRepClonotypes")
    repSeqObj = addItemToObject(repSeqObj, DAClonotypeAbundanceMatrix, 
        "subRepClonotypeAbundanceMatrix")
    return(repSeqObj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
