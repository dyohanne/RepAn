\name{findDAClusters}
\alias{findDAClusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findDAClusters(repSeqObj, abundanceType = c("cAbundance", "cRelAbundance", "cRelCloneSize"), testType = c("t.test", "wilcox.test", "RankProd"), paired = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repSeqObj}{
%%     ~~Describe \code{repSeqObj} here~~
}
  \item{abundanceType}{
%%     ~~Describe \code{abundanceType} here~~
}
  \item{testType}{
%%     ~~Describe \code{testType} here~~
}
  \item{paired}{
%%     ~~Describe \code{paired} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (repSeqObj, abundanceType = c("cAbundance", "cRelAbundance", 
    "cRelCloneSize"), testType = c("t.test", "wilcox.test", "RankProd"), 
    paired = F, ...) 
{
    grpLevels = levels(factor(repSeqObj$group))
    numberOfGroups = length(grpLevels)
    repSeqObj = getClusterAbundancesTable(repSeqObj)
    if (!abundanceType \%in\% c("cAbundance", "cRelAbundance", 
        "cRelCloneSize")) 
        stop("Abundance type for differential abundance testing not given. Please enter one of cAbundance, cRelAbundance or cRelCloneSize")
    if (abundanceType == "cAbundance") {
        selectedAbundanceTable = "cAbundanceTable"
    }
    if (abundanceType == "cRelAbundance") {
        selectedAbundanceTable = "cRelativeAbundanceTable"
    }
    if (abundanceType == "cRelCloneSize") {
        selectedAbundanceTable = "cRelativeAverageCloneSizeTable"
    }
    if (paired == T) {
        tooFewSubReps <- apply(repSeqObj[[selectedAbundanceTable]], 
            1, function(x) sum(!is.na(x[repSeqObj$group == grpLevels[2]]/x[repSeqObj$group == 
                grpLevels[1]])) > 2)
    }
    else {
        tooFewSubReps <- apply(repSeqObj[[selectedAbundanceTable]], 
            1, function(x) sum(!is.na(x[repSeqObj$group == grpLevels[1]])) > 
                2 & sum(!is.na(x[repSeqObj$group == grpLevels[2]])) > 
                2)
    }
    selectedSubrepertoiresTable <- repSeqObj[[selectedAbundanceTable]][tooFewSubReps, 
        ]
    repSeqObj = addItemToObject(repSeqObj, selectedSubrepertoiresTable, 
        "cSelectedSubRepertoireTable")
    selectedAbundanceTable = "cSelectedSubRepertoireTable"
    if (!testType \%in\% c("t.test", "wilcox.test", "RankProd")) {
        cat("Warning: Proper test type not given. RankProd test will be used.")
        testType = "RankProd"
    }
    else {
        if (testType[1] == "RankProd") {
            if (paired) {
                d = repSeqObj[[selectedAbundanceTable]][, repSeqObj$group == 
                  grpLevels[1]]/repSeqObj[[selectedAbundanceTable]][, 
                  repSeqObj$group == grpLevels[2]]
                tocl = rep(1, ncol(d))
            }
            else {
                d = repSeqObj[[selectedAbundanceTable]]
                tocl = c(rep(0, sum(repSeqObj$group == grpLevels[1])), 
                  rep(1, sum(repSeqObj$group == grpLevels[2])))
            }
            d_n = as.data.frame(normalize.quantiles(as.matrix(d)))
            colnames(d_n) <- colnames(d)
            rownames(d_n) <- rownames(d)
            d <- d_n
            sink("rpout")
            Rp.out <- RP(d, cl = tocl, logged = F, plot = F, 
                rand = 123, na.rm = T)
            sink(NULL)
            topSubReps = topGene(Rp.out, gene.names = rownames(d), 
                cutoff = 1, method = "pval", logged = F)
            if (nrow(topSubReps$Table1) > 0) {
                midx <- match(rownames(repSeqObj[[selectedAbundanceTable]]), 
                  rownames(topSubReps$Table1))
                midx <- midx[!is.na(midx)]
                repSeqObj[[selectedAbundanceTable]] <- repSeqObj[[selectedAbundanceTable]][which((rownames(repSeqObj[[selectedAbundanceTable]]) \%in\% 
                  rownames(topSubReps$Table1))), ]
                pval = topSubReps$Table1[midx, colnames(topSubReps$Table1) == 
                  "P.value"]
            }
            else {
            }
        }
        else {
            if (testType[1] == "t.test") {
                selectedTest = t.test
            }
            if (testType[1] == "wilcox.test") {
                selectedTest = wilcox.test
            }
            pval = apply(repSeqObj[[selectedAbundanceTable]], 
                1, function(x) selectedTest(x[repSeqObj$group == 
                  grpLevels[1]], x[repSeqObj$group == grpLevels[2]], 
                  paired = paired, na.rm = T, ...)$p.value)
        }
    }
    padj = round(p.adjust(pval, "BH"), 3)
    ratioAbundanceMeans = log2(rowMeans(repSeqObj[[selectedAbundanceTable]][, 
        repSeqObj$group == grpLevels[2]], na.rm = T)/rowMeans(repSeqObj[[selectedAbundanceTable]][, 
        repSeqObj$group == grpLevels[1]], na.rm = T))
    if (testType[1] == "RankProd") 
        ratioAbundanceMeans = log2(1/topSubReps$Table1[rownames(repSeqObj[[selectedAbundanceTable]]), 
            colnames(topSubReps$Table1) == "FC:(class1/class2)"])
    cDaResult = data.frame(logFC = ratioAbundanceMeans, pval, 
        padj)
    repSeqObj = addItemToObject(repSeqObj, cDaResult, "cDaResult")
    return(repSeqObj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
