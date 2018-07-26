\name{compareAbundanceInPairedSamples}
\alias{compareAbundanceInPairedSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compareAbundanceInPairedSamples(freqTable, pcutoff = 0.01, pairs, writeResult = F, resultDir = "matcheDA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{freqTable}{
%%     ~~Describe \code{freqTable} here~~
}
  \item{pcutoff}{
%%     ~~Describe \code{pcutoff} here~~
}
  \item{pairs}{
%%     ~~Describe \code{pairs} here~~
}
  \item{writeResult}{
%%     ~~Describe \code{writeResult} here~~
}
  \item{resultDir}{
%%     ~~Describe \code{resultDir} here~~
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
function (freqTable, pcutoff = 0.01, pairs, writeResult = F, 
    resultDir = "matcheDA") 
{
    sam = freqTable
    resultlist <- list()
    pvals = c()
    ORs = c()
    CIs = c()
    print(paste("There are", ncol(sam), "samples in the data. Comparison is done for pair of samples..."))
    print(paste("Number of features to compare : ", nrow(sam)))
    dir.create(file.path(resultDir), showWarnings = FALSE)
    rnds = length(pairs)/2
    for (k in 1:rnds) {
        kpair = k + rnds
        res = compareAbundanceInTwoSamplesForShared(sam[, c(k, 
            kpair)])
        colnames(res) <- c(paste(colnames(sam[, c(k, kpair)])[1], 
            colnames(sam[, c(k, kpair)])[2], "padjusted", sep = "_"), 
            paste(colnames(sam[, c(k, kpair)])[1], colnames(sam[, 
                c(k, kpair)])[2], "OR", sep = "_"), paste(colnames(sam[, 
                c(k, kpair)])[1], colnames(sam[, c(k, kpair)])[2], 
                "CI", sep = "_"))
        print(str(res[, 1]))
        siglist = as.numeric(res[, 1]) < pcutoff
        resSig = res[siglist, ]
        print(head(siglist))
        resultlist[[colnames(sam[, c(k, kpair)])[1]]] <- cbind(sam[rownames(resSig), 
            c(k, kpair)], resSig)
        if (writeResult == T) {
            print("Result written to file....")
            write.table(res, file = paste(resultDir, "/", colnames(sam[, 
                c(k, kpair)])[1], "_", colnames(sam[, c(k, kpair)])[2], 
                "_Abundance_fishers_comparison_result_forPairSamples.txt", 
                sep = ""), row.names = T, sep = "\t")
        }
    }
    return(resultlist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
