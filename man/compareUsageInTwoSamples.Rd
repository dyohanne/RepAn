\name{compareUsageInTwoSamples}
\alias{compareUsageInTwoSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compareUsageInTwoSamples(freqTable, pcutoff = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{freqTable}{
%%     ~~Describe \code{freqTable} here~~
}
  \item{pcutoff}{
%%     ~~Describe \code{pcutoff} here~~
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
function (freqTable, pcutoff = 0.05) 
{
    sam = freqTable
    pvals = c()
    ORs = c()
    CIs = c()
    print(head(sam))
    for (i in 1:nrow(sam)) {
        current <- sam[i, ]
        othersSum <- colSums(sam[-i, ])
        dataForComparison = matrix(c(as.numeric(current), as.numeric(othersSum)), 
            nrow = 2, byrow = T)
        fresult = fisher.test(dataForComparison)
        pvals = c(pvals, fresult$p.value)
        ORs = c(ORs, as.numeric(fresult$estimate))
        CIs = c(CIs, paste(fresult$conf[1], fresult$conf[2], 
            sep = "-"))
    }
    pvals = p.adjust(pvals, "BH")
    res = cbind(pvals, ORs, CIs)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
