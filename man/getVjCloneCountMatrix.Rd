\name{getVjCloneCountMatrix}
\alias{getVjCloneCountMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getVjCloneCountMatrix(reps, vj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reps}{
%%     ~~Describe \code{reps} here~~
}
  \item{vj}{
%%     ~~Describe \code{vj} here~~
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
function (reps, vj) 
{
    m <- vj
    colnames(m) <- c("VGENENAME", "JGENENAME")
    for (rep in reps) {
        print(rep)
        p <- get(rep)
        p <- getProductive(p)
        vjPairAgg <- aggregate(AMINOACID ~ VGENENAME + JGENENAME, 
            data = p, FUN = length)
        print(head(vjPairAgg))
        colnames(vjPairAgg) <- c("VGENENAME", "JGENENAME", rep)
        m <- merge(m, vjPairAgg, by = c("VGENENAME", "JGENENAME"), 
            all.x = T)
    }
    m[is.na(m)] <- 0
    m[, 1] <- paste(m[, 1], m[, 2])
    m <- m[, -2]
    rownames(m) <- m[, 1]
    m <- m[, -1]
    return(m)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
