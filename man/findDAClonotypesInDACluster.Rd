\name{findDAClonotypesInDACluster}
\alias{findDAClonotypesInDACluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findDAClonotypesInDACluster(clonotypes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clonotypes}{
%%     ~~Describe \code{clonotypes} here~~
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
function (clonotypes) 
{
    clonotypes = as.vector(clonotypes)
    imp_clones_qiao <- read.table("Important clones qiao.txt", 
        header = T, sep = "\t", dec = ",")
    imp_clones_qiao = unique(as.vector(imp_clones_qiao[, 1]))
    imp_clones_arnold1 <- read.table("Important clones Arnold_CD4_gluten reactive.txt", 
        header = T, sep = "\t", dec = ",")
    imp_clones_arnold1 = unique(as.vector(imp_clones_arnold1[, 
        1]))
    imp_clones_arnold2 <- read.table("Important clones Arnold_CD8_isolated from challenged patients.txt", 
        header = T, sep = "\t", dec = ",")
    imp_clones_arnold2 = unique(as.vector(imp_clones_arnold2[, 
        1]))
    imp_clones_arnold_total = c(as.vector(imp_clones_arnold1), 
        as.vector(imp_clones_arnold2))
    imp_clones_petersen <- read.table("Important clones Petersen.txt", 
        header = T, sep = "\t", dec = ",")
    imp_clones_petersen = unique(as.vector(imp_clones_petersen[, 
        8]))
    knownClones <- unique(c(imp_clones_qiao, imp_clones_arnold1, 
        imp_clones_arnold2, imp_clones_petersen))
    found_in_knownclones = clonotypes[clonotypes \%in\% knownClones]
    return(found_in_knownclones)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
