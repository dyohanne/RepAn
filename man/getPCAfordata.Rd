\name{getPCAfordata}
\alias{getPCAfordata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPCAfordata(dd, groups, scaling = T, oprefix = "O")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dd}{
%%     ~~Describe \code{dd} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
}
  \item{scaling}{
%%     ~~Describe \code{scaling} here~~
}
  \item{oprefix}{
%%     ~~Describe \code{oprefix} here~~
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
function (dd, groups, scaling = T, oprefix = "O") 
{
    dd = t(dd)
    pr <- prcomp(dd, scale = scaling)
    print(summary(pr))
    sampleNames = rownames(dd)
    s = summary(pr)$importance[2, ]
    flname = paste(oprefix, "PCA.pdf")
    pdf(file = flname, useDingbats = FALSE)
    tobedrawn = data.frame(pr$x, group = factor(groups), snames = rownames(dd))
    write.table(pr$x, file = paste(oprefix, "_PCArotatedData.txt"), 
        row.names = T, col.names = T, sep = "\t")
    pcaP = ggplot(tobedrawn, aes(x = PC1, y = PC2, shape = group, 
        colour = snames)) + xlab(paste("PC1 (", format(round(s[1] * 
        100), nsmall = 2), "\%)")) + ylab(paste("PC2 (", format(round(s[2] * 
        100), nsmall = 2), "\%)")) + geom_point(size = 6) + scale_fill_brewer(palette = "Set1")
    pcaP = pcaP + theme_bw() + theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank()) + theme(axis.line = element_line(colour = "black"))
    print(pcaP)
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
