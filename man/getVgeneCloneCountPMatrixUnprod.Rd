\name{getVgeneCloneCountPMatrixUnprod}
\alias{getVgeneCloneCountPMatrixUnprod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getVgeneCloneCountPMatrixUnprod(reps, g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reps}{
%%     ~~Describe \code{reps} here~~
}
  \item{g}{
%%     ~~Describe \code{g} here~~
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
function (reps, g) 
{
    m <- as.data.frame(g)
    colnames(m) <- c("VGENENAME")
    for (rep in reps) {
        print(rep)
        p <- get(rep)
        p <- getUnproductive(p)
        print(nrow(p))
        vfamAgg <- aggregate(AMINOACID ~ VGENENAME, data = p, 
            FUN = length)
        vfamAgg[, 2] <- vfamAgg[, 2]/sum(vfamAgg[, 2])
        colnames(vfamAgg) <- c("VGENENAME", paste(rep))
        m <- merge(m, vfamAgg, by = "VGENENAME", all.x = T)
    }
    m[is.na(m)] <- 0
    rownames(m) <- m[, 1]
    m <- m[, -1]
    return(m)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
