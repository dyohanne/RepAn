\name{compareUsageProportionsPaired}
\alias{compareUsageProportionsPaired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compareUsageProportionsPaired(freqTable, pcutoff = 0.05, pairs)
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
function (freqTable, pcutoff = 0.05, pairs) 
{
    sam = freqTable
    resultTable <- sam
    pvals = c()
    ORs = c()
    CIs = c()
    print(paste("There are", ncol(sam), "samples in the data. Comparison is done for pair of samples..."))
    print(paste("Number of features to compare : ", nrow(sam)))
    rnds = length(pairs)/2
    for (k in 1:rnds) {
        kpair = k + rnds
        res = compareUsageInTwoSamples(sam[, c(k, kpair)])
        colnames(res) <- c(paste(colnames(sam[, c(k, kpair)])[1], 
            colnames(sam[, c(k, kpair)])[2], "padjusted", sep = "_"), 
            paste(colnames(sam[, c(k, kpair)])[1], colnames(sam[, 
                c(k, kpair)])[2], "OR", sep = "_"), paste(colnames(sam[, 
                c(k, kpair)])[1], colnames(sam[, c(k, kpair)])[2], 
                "CI", sep = "_"))
        resultTable <- cbind(resultTable, res)
    }
    print("Result written to file....")
    write.table(resultTable, file = "Usage_fishers_comparison_result_forPairSamples.txt", 
        row.names = T, sep = "\t")
    return(resultTable)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
