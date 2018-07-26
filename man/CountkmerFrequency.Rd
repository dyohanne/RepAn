\name{CountkmerFrequency}
\alias{CountkmerFrequency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CountkmerFrequency(seqs, type = "NT", k = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqs}{
%%     ~~Describe \code{seqs} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
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
function (seqs, type = "NT", k = 4) 
{
    DNAletters = c("A", "C", "G", "T")
    AAletters = c("A", "C", "D", "E", "F", "G", "H", "I", "K", 
        "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
    if (type == "NT") {
        allCombinations = do.call(expand.grid, rep(list(DNAletters), 
            k))
        kmers = apply(allCombinations, 1, paste, collapse = "")
    }
    else {
        allCombinations = do.call(expand.grid, rep(list(AAletters), 
            k))
        kmers = apply(allCombinations, 1, paste, collapse = "")
    }
    kmerCount <- sapply(seqs, function(x) {
        sapply(kmers, function(y) {
            sum(gregexpr(paste("(?=", y, ")", sep = ""), x, perl = T)[[1]] > 
                0)
        })
    })
    kmerCount <- t(kmerCount)
    kmerCount <- kmerCount[, order(colnames(kmerCount))]
    return(kmerCount)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
