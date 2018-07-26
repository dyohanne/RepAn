\name{getClonalKmerNtFrequencies}
\alias{getClonalKmerNtFrequencies}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getClonalKmerNtFrequencies(reps, selectionType = "random", nClones = 100, kMer = 4, asProb = F, vGenes = NULL, seqStatus = "prod")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reps}{
%%     ~~Describe \code{reps} here~~
}
  \item{selectionType}{
%%     ~~Describe \code{selectionType} here~~
}
  \item{nClones}{
%%     ~~Describe \code{nClones} here~~
}
  \item{kMer}{
%%     ~~Describe \code{kMer} here~~
}
  \item{asProb}{
%%     ~~Describe \code{asProb} here~~
}
  \item{vGenes}{
%%     ~~Describe \code{vGenes} here~~
}
  \item{seqStatus}{
%%     ~~Describe \code{seqStatus} here~~
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
function (reps, selectionType = "random", nClones = 100, kMer = 4, 
    asProb = F, vGenes = NULL, seqStatus = "prod") 
{
    kmerFreqTable = NULL
    repnames = c()
    for (rep in reps) {
        print(rep)
        p <- get(rep)
        if (seqStatus == "prod") {
            p <- getProductive(p)
        }
        else {
            p <- getUnproductive(p)
        }
        if (!is.null(vGenes)) {
            p <- p[p$VGENENAME \%in\% vGenes, ]
        }
        if (selectionType == "random") {
            sampling.index = sample(1:nrow(p), nClones, replace = F)
        }
        else {
            sampling.index = 1:nClones
        }
        p <- p[sampling.index, ]
        print(nrow(p))
        seqs <- unique(p$NUCLEOTIDE)
        seqSet <- DNAStringSet(seqs)
        seq_mers <- oligonucleotideFrequency(seqSet, width = kMer, 
            as.prob = asProb)
        kmerFreqTable = rbind(kmerFreqTable, seq_mers)
        repnames = c(repnames, rep(rep, nClones))
    }
    rownames(kmerFreqTable) <- repnames
    kmerFreqTable <- t(kmerFreqTable)
    print("Writing result to file....")
    write.table(kmerFreqTable, file = "clonal4merFrequencyTable.txt", 
        row.names = T, col.names = T, sep = "\t")
    return(list(freqTable = kmerFreqTable, sams = repnames))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
