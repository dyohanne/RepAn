\name{readSample}
\alias{readSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
readSample(samplename, cloneStatus = "p")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samplename}{
%%     ~~Describe \code{samplename} here~~
}
  \item{cloneStatus}{
%%     ~~Describe \code{cloneStatus} here~~
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
function (samplename, cloneStatus = "p") 
{
    sam <- read.table(samplename, header = T, sep = "\t", dec = ".", 
        stringsAsFactors = F)
    colnames(sam) <- toupper(colnames(sam))
    sam <- sam[, which(colnames(sam) == "NUCLEOTIDE"):length(colnames(sam))]
    CDR3NT = substring(sam$NUCLEOTIDE, sam$VINDEX + 1, sam$VINDEX + 
        sam$CDR3LENGTH)
    sam <- data.frame(sam, CDR3NT, stringsAsFactors = F)
    if (sum(colnames(sam) == "COPY") == 1) {
        colnames(sam)[which(colnames(sam) == "COPY")] <- "UNNORMALIZEDCOUNT"
        colnames(sam)[which(colnames(sam) == "NORMALIZEDCOPY")] <- "COUNT"
        colnames(sam)[which(colnames(sam) == "NORMALIZEDFREQUENCY")] <- "FREQUENCYCOUNT"
    }
    if (cloneStatus == "p") {
        sam = getProductive(sam)
    }
    else if (cloneStatus == "unp") {
        sam = getUnproductive(sam)
    }
    return(sam)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
