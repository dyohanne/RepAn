\name{compareAbundanceInPairedSamplesForRanking}
\alias{compareAbundanceInPairedSamplesForRanking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compareAbundanceInPairedSamplesForRanking(samObj, freqTable, pairs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samObj}{
%%     ~~Describe \code{samObj} here~~
}
  \item{freqTable}{
%%     ~~Describe \code{freqTable} here~~
}
  \item{pairs}{
%%     ~~Describe \code{pairs} here~~
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
function (samObj, freqTable, pairs) 
{
    sam = round(freqTable)
    resultlist <- list()
    rnds = length(pairs)/2
    for (k in 1:rnds) {
        kpair = k + rnds
        pvals = c()
        ORs = c()
        CIs = c()
        inSample = c()
        NttoAA = c()
        for (i in 1:nrow(sam)) {
            current <- c()
            othersSum <- c()
            nttoaa1 <- length(samObj$sampleData[[k]][samObj$sampleData[[k]]$AMINOACID \%in\% 
                rownames(sam[i, , drop = F]), c("COUNT")])
            nttoaa2 <- length(samObj$sampleData[[kpair]][samObj$sampleData[[kpair]]$AMINOACID \%in\% 
                rownames(sam[i, , drop = F]), c("COUNT")])
            if (nttoaa1 == 0) 
                nttoaa1 = 1
            NttoAA <- c(NttoAA, nttoaa2/nttoaa1)
            current[1] <- sum(samObj$sampleData[[k]][samObj$sampleData[[k]]$AMINOACID \%in\% 
                rownames(sam[i, , drop = F]), c("COUNT")])
            current[2] <- sum(samObj$sampleData[[kpair]][samObj$sampleData[[kpair]]$AMINOACID \%in\% 
                rownames(sam[i, , drop = F]), c("COUNT")])
            current[is.na(current)] <- 0
            if (current[1] == 0 & current[2] == 0) {
                inSample <- c(inSample, F)
                othersSum[1] <- sum(samObj$sampleData[[k]][!samObj$sampleData[[k]]$AMINOACID \%in\% 
                  rownames(sam[i, , drop = F]), c("COUNT")])
                othersSum[2] <- sum(samObj$sampleData[[kpair]][!samObj$sampleData[[kpair]]$AMINOACID \%in\% 
                  rownames(sam[i, , drop = F]), c("COUNT")])
            }
            else if (current[1] == 0 & current[2] != 0) {
                inSample <- c(inSample, T)
                othersSum[1] <- sum(samObj$sampleData[[k]][, 
                  c("COUNT")])
                othersSum[2] <- sum(samObj$sampleData[[kpair]][!samObj$sampleData[[kpair]]$AMINOACID \%in\% 
                  rownames(sam[i, , drop = F]), c("COUNT")])
            }
            else if (current[1] != 0 & current[2] == 0) {
                inSample <- c(inSample, T)
                othersSum[1] <- sum(samObj$sampleData[[k]][!samObj$sampleData[[k]]$AMINOACID \%in\% 
                  rownames(sam[i, , drop = F]), c("COUNT")])
                othersSum[2] <- sum(samObj$sampleData[[kpair]][, 
                  c("COUNT")])
            }
            else {
                inSample <- c(inSample, T)
                othersSum[1] <- sum(samObj$sampleData[[k]][!samObj$sampleData[[k]]$AMINOACID \%in\% 
                  rownames(sam[i, , drop = F]), c("COUNT")])
                othersSum[2] <- sum(samObj$sampleData[[kpair]][!samObj$sampleData[[kpair]]$AMINOACID \%in\% 
                  rownames(sam[i, , drop = F]), c("COUNT")])
            }
            dataForComparison = matrix(c(as.numeric(current), 
                as.numeric(othersSum)), nrow = 2, byrow = T)
            dataForComparison[is.na(dataForComparison)] <- 0
            dataForComparison = dataForComparison + 1
            dataForComparison <- round(dataForComparison)
            fresult = fisher.test(dataForComparison[, c(2, 1)])
            pvals = c(pvals, fresult$p.value)
            ORs = c(ORs, as.numeric(fresult$estimate))
            CIs = c(CIs, paste(fresult$conf[1], fresult$conf[2], 
                sep = "-"))
        }
        res = cbind(inSample * 1, NttoAA, as.numeric(pvals), 
            ORs, CIs)
        colnames(res) <- c("inSample", "NTtoAA", paste(colnames(sam[, 
            c(k, kpair)])[1], colnames(sam[, c(k, kpair)])[2], 
            "pvalue", sep = "_"), paste(colnames(sam[, c(k, kpair)])[1], 
            colnames(sam[, c(k, kpair)])[2], "OR", sep = "_"), 
            paste(colnames(sam[, c(k, kpair)])[1], colnames(sam[, 
                c(k, kpair)])[2], "CI", sep = "_"))
        resultlist[[k]] <- res
    }
    avPvals <- c()
    avORs <- c()
    avNTtoAAs <- c()
    nOfSampesForClones <- c()
    for (i in 1:nrow(sam)) {
        nSamplesWithTheClone <- as.numeric(sapply(resultlist, 
            function(x) x[i, 1]))
        NttoasInsamples <- as.numeric(sapply(resultlist, function(x) x[i, 
            2]))
        pvals <- as.numeric(sapply(resultlist, function(x) x[i, 
            3]))
        ors <- as.numeric(sapply(resultlist, function(x) x[i, 
            4]))
        mpval <- mean(pvals[nSamplesWithTheClone == 1])
        mors <- mean(ors[nSamplesWithTheClone == 1])
        mNTtoAA <- mean(NttoasInsamples[nSamplesWithTheClone == 
            1])
        avPvals <- c(avPvals, mpval)
        avORs <- c(avORs, mors)
        avNTtoAAs <- c(avNTtoAAs, mNTtoAA)
        nOfSampesForClones <- c(nOfSampesForClones, sum(sam[i, 
            pairs == 1] > 0) - sum(sam[i, pairs == 0] > 0))
    }
    toRet <- data.frame(avNTtoAAs, avPvals, avORs, nOfSampesForClones)
    rownames(toRet) <- rownames(sam)
    return(toRet)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
