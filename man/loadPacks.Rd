\name{loadPacks}
\alias{loadPacks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loadPacks(package.list = c("ggplot2", "seqRFLP", "stringr", "permute", "cluster", "data.table", "NbClust", "doParallel", "dendextend", "kmer", "ape", "Matrix", "wordspace", "dynamicTreeCut", "e1071", "fclust", "randomForest", "preprocessCore", "gplots"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{package.list}{
%%     ~~Describe \code{package.list} here~~
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
function (package.list = c("ggplot2", "seqRFLP", "stringr", "permute", 
    "cluster", "data.table", "NbClust", "doParallel", "dendextend", 
    "kmer", "ape", "Matrix", "wordspace", "dynamicTreeCut", "e1071", 
    "fclust", "randomForest", "preprocessCore", "gplots")) 
{
    new.packages <- package.list[!(package.list \%in\% installed.packages()[, 
        "Package"])]
    if (length(new.packages)) 
        try(install.packages(new.packages))
    lapply(eval(package.list), require, character.only = TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
