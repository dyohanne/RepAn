\name{determineWeight}
\alias{determineWeight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
determineWeight(kmer, seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kmer}{
%%     ~~Describe \code{kmer} here~~
}
  \item{seq}{
%%     ~~Describe \code{seq} here~~
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
function (kmer, seq) 
{
    kmerPositions = as.numeric(gregexpr(kmer, seq)[[1]])
    seqLength <- nchar(seq)
    weightGroups <- seq(seqLength/3, seqLength, seqLength/3)
    kmerWts = c()
    for (i in kmerPositions) {
        if (i == -1) {
            wt = 0
        }
        else {
            wps = sum(weightGroups > i)
            if (wps == 3) {
                wt = 5
            }
            if (wps == 2) {
                wt = 10
            }
            if (wps == 1) {
                wt = 1
            }
        }
        kmerWts = c(kmerWts, wt)
    }
    return(sum(kmerWts))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
