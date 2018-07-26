\name{findOptimalClusters}
\alias{findOptimalClusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findOptimalClusters(repSeqObj, k, clusterby = "NT", kmerWidth = 4, posWt = F, distMethod = "euclidean", useDynamicTreeCut = T, matchingMethod = c("hc", "km", "og"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repSeqObj}{
%%     ~~Describe \code{repSeqObj} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{clusterby}{
%%     ~~Describe \code{clusterby} here~~
}
  \item{kmerWidth}{
%%     ~~Describe \code{kmerWidth} here~~
}
  \item{posWt}{
%%     ~~Describe \code{posWt} here~~
}
  \item{distMethod}{
%%     ~~Describe \code{distMethod} here~~
}
  \item{useDynamicTreeCut}{
%%     ~~Describe \code{useDynamicTreeCut} here~~
}
  \item{matchingMethod}{
%%     ~~Describe \code{matchingMethod} here~~
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
function (repSeqObj, k, clusterby = "NT", kmerWidth = 4, posWt = F, 
    distMethod = "euclidean", useDynamicTreeCut = T, matchingMethod = c("hc", 
        "km", "og")) 
{
    withinSampleClusters = list()
    for (i in 1:length(repSeqObj$samNames)) {
        cResult = getClusterLables(repSeqObj$scaledSampleData[[i]], 
            k, clusterby, kmerWidth, posWt, distMethod, useDynamicTreeCut = T)
        withinSampleClusters[[length(withinSampleClusters) + 
            1]] <- cResult
    }
    names(withinSampleClusters) <- repSeqObj$samNames
    repSeqObj <- addItemToObject(repSeqObj, withinSampleClusters, 
        "withinSampleClusters")
    if (matchingMethod == "og") {
        matchTable <- getClusterMatches(repSeqObj, matchingMethod = "og")
    }
    else if (matchingMethod == "km") {
        matchTable <- getClusterMatches(repSeqObj, matchingMethod = "km")
    }
    else if (matchingMethod == "hc") {
        matchTable <- getClusterMatches(repSeqObj, matchingMethod = "hc", 
            distMethod = distMethod)
    }
    else {
        cat("Matching method can only be one of hc, km or og.hc method will be used.\n")
        matchTable <- getClusterMatches(repSeqObj, matchingMethod = "hc", 
            distMethod = distMethod)
    }
    repSeqObj <- addItemToObject(repSeqObj, matchTable, "clusterMatchTable")
    return(repSeqObj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
