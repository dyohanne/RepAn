\name{compareClusterAbundancesPaired}
\alias{compareClusterAbundancesPaired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compareClusterAbundancesPaired(sam1, sam2, clusMatchTable, s1cls, s2cls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sam1}{
%%     ~~Describe \code{sam1} here~~
}
  \item{sam2}{
%%     ~~Describe \code{sam2} here~~
}
  \item{clusMatchTable}{
%%     ~~Describe \code{clusMatchTable} here~~
}
  \item{s1cls}{
%%     ~~Describe \code{s1cls} here~~
}
  \item{s2cls}{
%%     ~~Describe \code{s2cls} here~~
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
function (sam1, sam2, clusMatchTable, s1cls, s2cls) 
{
    clusterFcs = c()
    clsFCsObservedPermuted = NULL
    for (ck in 1:nrow(clusMatchTable)) {
        consensusK <- clusMatchTable[ck, ]
        sam1selected <- sam1[s1cls$clusters == consensusK[1], 
            ]
        sam2selected <- sam2[s2cls$clusters == consensusK[2], 
            ]
        sam1CountPerClone <- sum(sam1$COUNT)/nrow(sam1)
        sam2CountPerClone <- sum(sam2$COUNT)/nrow(sam2)
        sam1SelectedCountPerClone <- sum(sam1selected$COUNT)/nrow(sam1selected)
        sam2SelectedCountPerClone <- sum(sam2selected$COUNT)/nrow(sam2selected)
        sam1SelectedVsTotal <- sam1SelectedCountPerClone/sam1CountPerClone
        sam2SelectedVsTotal <- sam2SelectedCountPerClone/sam2CountPerClone
        sam1vssam2Selected <- round(sam2SelectedVsTotal/sam1SelectedVsTotal, 
            digits = 2)
        clusterFcs <- c(clusterFcs, sam1vssam2Selected)
    }
    sam1new = sam1cls
    sam2new = sam2cls
    clsFCsObservedPermuted = rbind(clsFCsObservedPermuted, clusterFcs)
    for (s in 1:100) {
        sam1new$clusters <- sample(as.numeric(sam1cls$clusters), 
            length(sam1cls$clusters), replace = F)
        sam2new$clusters <- sample(as.numeric(sam2cls$clusters), 
            length(sam2cls$clusters), replace = F)
        pClusterFc <- getClusterFoldChanges(sam1, sam2, clusMatchTable, 
            sam1new, sam2new)
        clsFCsObservedPermuted <- rbind(clsFCsObservedPermuted, 
            pClusterFc)
    }
    ps <- c()
    for (i in 1:ncol(clsFCsObservedPermuted)) {
        obs <- clsFCsObservedPermuted[1, i]
        otherValues <- clsFCsObservedPermuted[-1, i]
        pval <- sum(abs(otherValues) >= abs(obs))/nrow(clsFCsObservedPermuted)
        ps <- c(ps, pval)
    }
    return(list(observedClusterFC = clusterFcs, Pval = ps))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
