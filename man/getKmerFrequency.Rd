\name{getKmerFrequency}
\alias{getKmerFrequency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getKmerFrequency(seqs, type = "NT", k = 4, normForLength = F)
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
  \item{normForLength}{
%%     ~~Describe \code{normForLength} here~~
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
function (seqs, type = "NT", k = 4, normForLength = F) 
{
    if (type == "NT") {
        seqSet <- DNAStringSet(seqs)
        seq_mers <- oligonucleotideFrequency(seqSet, width = k, 
            as.prob = F)
        rownames(seq_mers) <- seqs
    }
    else {
        seqList <- as.list(seqs)
        seqAAbin <- as.AAbin(seqList)
        seq_mers <- kcount(seqAAbin, k = k)
        rownames(seq_mers) <- seqs
    }
    if (normForLength == T) 
        seq_mers <- t(apply(seq_mers, 1, function(x) x/sum(x)))
    seq_mers <- as(seq_mers, "sparseMatrix")
    return(seq_mers)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
