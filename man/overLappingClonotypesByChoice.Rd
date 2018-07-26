\name{overLappingClonotypesByChoice}
\alias{overLappingClonotypesByChoice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
overLappingClonotypesByChoice(reps, cutoff = 100, cutoffby = "perc", seqType = "AA", unProdStatus = "unprod")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reps}{
%%     ~~Describe \code{reps} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
}
  \item{cutoffby}{
%%     ~~Describe \code{cutoffby} here~~
}
  \item{seqType}{
%%     ~~Describe \code{seqType} here~~
}
  \item{unProdStatus}{
%%     ~~Describe \code{unProdStatus} here~~
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
function (reps, cutoff = 100, cutoffby = "perc", seqType = "AA", 
    unProdStatus = "unprod") 
{
    cutoffby = toupper(cutoffby)
    m <- as.data.frame(combinations(n = length(reps), r = 2, 
        v = reps, repeats.allowed = TRUE))
    colnames(m) <- c("Sample1", "Sample2")
    toadd = data.frame()
    unProdStatus <- toupper(unProdStatus)
    seqType <- toupper(seqType)
    if (seqType == "AA") {
        seqType <- "AMINOACID"
    }
    else {
        seqType <- "NUCLEOTIDE"
    }
    for (r in 1:nrow(m)) {
        overlapValuesProd = c()
        overlapValuesUnprod = c()
        overlapValuesTotalSam1Prod = c()
        overlapValuesTotalSam2Prod = c()
        overlapValuesTotalSam1Unprod = c()
        overlapValuesTotalSam2Unprod = c()
        overlapValuesProdJI = c()
        overlapValuesUnprodJI = c()
        overlapValuesProdJI_abundancebased = c()
        overlapValuesUnprodJI_abundancebased = c()
        overlapValuesProdJI_abundancebased2 = c()
        overlapValuesUnprodJI_abundancebased2 = c()
        NinsertionProd1 = c()
        NinsertionUnprod1 = c()
        NinsertionProd2 = c()
        NinsertionUnprod2 = c()
        print(r)
        sam1 = m[r, 1]
        sam2 = m[r, 2]
        samplesCompared = paste(sam1, sam2, sep = "_")
        sam1 <- get(as.character(sam1))
        sam2 <- get(as.character(sam2))
        sam1Prod <- getProductive(sam1)
        if (unProdStatus == "STOP") {
            sam1Unprod <- getHasStop(sam1)
        }
        else if (unProdStatus == "OUT") {
            sam1Unprod <- getOutOfFrame(sam1)
        }
        else {
            sam1Unprod <- getUnproductive(sam1)
        }
        sam2Prod <- getProductive(sam2)
        if (unProdStatus == "STOP") {
            sam2Unprod <- getHasStop(sam2)
        }
        else if (unProdStatus == "OUT") {
            sam2Unprod <- getOutOfFrame(sam2)
        }
        else {
            sam2Unprod <- getUnproductive(sam2)
        }
        indSam1 <- which(colnames(sam1) == seqType)
        indSam2 <- which(colnames(sam2) == seqType)
        if (cutoffby == "PERC") {
            topclones_sam1Prod = getAbundantClonotypesbyperc(sam1Prod, 
                cutoff)
            topclones_sam1Unprod = getAbundantClonotypesbyperc(sam1Unprod, 
                cutoff)
            topclones_sam2Prod = getAbundantClonotypesbyperc(sam2Prod, 
                cutoff)
            topclones_sam2Unprod = getAbundantClonotypesbyperc(sam2Unprod, 
                cutoff)
        }
        else {
            topclones_sam1Prod = getAbundantClonotypesbynumber(sam1Prod, 
                cutoff)
            topclones_sam1Unprod = getAbundantClonotypesbynumber(sam1Unprod, 
                cutoff)
            topclones_sam2Prod = getAbundantClonotypesbynumber(sam2Prod, 
                cutoff)
            topclones_sam2Unprod = getAbundantClonotypesbynumber(sam2Unprod, 
                cutoff)
        }
        commonClonesProd = length(intersect(topclones_sam1Prod[, 
            indSam1], topclones_sam2Prod[, indSam2]))
        unionProd = length(union(topclones_sam1Prod[, indSam1], 
            topclones_sam2Prod[, indSam2]))
        commonClonesUnprod = length(intersect(topclones_sam1Unprod[, 
            indSam1], topclones_sam2Unprod[, indSam2]))
        unionUnprod = length(union(topclones_sam1Unprod[, indSam1], 
            topclones_sam2Unprod[, indSam2]))
        if (unionUnprod == 0) {
            unionUnprod = 1
        }
        if (unionProd == 0) {
            unionProd = 1
        }
        overlapValuesProd = c(overlapValuesProd, as.numeric(commonClonesProd))
        overlapValuesUnprod = c(overlapValuesUnprod, as.numeric(commonClonesUnprod))
        overlapValuesProdJI = c(overlapValuesProdJI, commonClonesProd/unionProd)
        overlapValuesUnprodJI = c(overlapValuesUnprodJI, commonClonesUnprod/unionUnprod)
        commonClonesProd_U_data = topclones_sam1Prod[topclones_sam1Prod[, 
            indSam1] \%in\% intersect(topclones_sam1Prod[, indSam1], 
            topclones_sam2Prod[, indSam2]), ]
        commonClonesProd_U = sum(commonClonesProd_U_data$COUNT)/sum(topclones_sam1Prod$COUNT)
        commonClonesProd_U_abundance = sum(commonClonesProd_U_data$COUNT)
        commonClonesProd_V_data = topclones_sam2Prod[topclones_sam2Prod[, 
            indSam2] \%in\% intersect(topclones_sam1Prod[, indSam1], 
            topclones_sam2Prod[, indSam2]), ]
        commonClonesProd_V = sum(commonClonesProd_V_data$COUNT)/sum(topclones_sam2Prod$COUNT)
        commonClonesProd_V_abundance = sum(commonClonesProd_V_data$COUNT)
        overlapValuesTotalSam1Prod = c(overlapValuesTotalSam1Prod, 
            as.numeric(commonClonesProd_U_abundance))
        overlapValuesTotalSam2Prod = c(overlapValuesTotalSam2Prod, 
            as.numeric(commonClonesProd_V_abundance))
        uv = commonClonesProd_U * commonClonesProd_V
        uplusv = commonClonesProd_U + commonClonesProd_V
        AbundanceBasedJIProd = uv/(uplusv - uv)
        overlapValuesProdJI_abundancebased = c(overlapValuesProdJI_abundancebased, 
            as.numeric(AbundanceBasedJIProd))
        totalSam1Sam2Prod = sum(topclones_sam1Prod$COUNT) + sum(topclones_sam2Prod$COUNT)
        AbundanceBasedJIProd2 = (commonClonesProd_U_abundance + 
            commonClonesProd_V_abundance)/totalSam1Sam2Prod
        overlapValuesProdJI_abundancebased2 = c(overlapValuesProdJI_abundancebased2, 
            as.numeric(AbundanceBasedJIProd2))
        commonClonesUnProd_U_data = topclones_sam1Unprod[topclones_sam1Unprod[, 
            indSam1] \%in\% intersect(topclones_sam1Unprod[, indSam1], 
            topclones_sam2Unprod[, indSam2]), ]
        commonClonesUnProd_U = sum(commonClonesUnProd_U_data$COUNT)/sum(topclones_sam1Unprod$COUNT)
        commonClonesUnProd_U_abundance = sum(commonClonesUnProd_U_data$COUNT)
        commonClonesUnProd_V_data = topclones_sam2Unprod[topclones_sam2Unprod[, 
            indSam2] \%in\% intersect(topclones_sam1Unprod[, indSam1], 
            topclones_sam2Unprod[, indSam2]), ]
        commonClonesUnProd_V = sum(commonClonesUnProd_V_data$COUNT)/sum(topclones_sam2Unprod$COUNT)
        commonClonesUnProd_V_abundance = sum(commonClonesUnProd_V_data$COUNT)
        overlapValuesTotalSam1Unprod = c(overlapValuesTotalSam1Unprod, 
            as.numeric(commonClonesUnProd_U_abundance))
        overlapValuesTotalSam2Unprod = c(overlapValuesTotalSam2Unprod, 
            as.numeric(commonClonesUnProd_V_abundance))
        uv_UnProd = commonClonesUnProd_U * commonClonesUnProd_V
        uplusv_UnProd = commonClonesUnProd_U + commonClonesUnProd_V
        AbundanceBasedJIUnProd = uv_UnProd/(uplusv_UnProd - uv_UnProd)
        overlapValuesUnprodJI_abundancebased = c(overlapValuesUnprodJI_abundancebased, 
            as.numeric(AbundanceBasedJIUnProd))
        totalSam1Sam2Unprod = sum(topclones_sam1Unprod$COUNT) + 
            sum(topclones_sam2Unprod$COUNT)
        AbundanceBasedJIUnProd2 = (commonClonesUnProd_U_abundance + 
            commonClonesUnProd_V_abundance)/totalSam1Sam2Unprod
        overlapValuesUnprodJI_abundancebased2 = c(overlapValuesUnprodJI_abundancebased2, 
            as.numeric(AbundanceBasedJIUnProd2))
        commonClonotypes1 = topclones_sam1Prod[topclones_sam1Prod[, 
            indSam1] \%in\% intersect(topclones_sam1Prod[, indSam1], 
            topclones_sam2Prod[, indSam1]), ]
        commonClonotypes2 = topclones_sam2Prod[topclones_sam2Prod[, 
            indSam2] \%in\% intersect(topclones_sam1Prod[, indSam1], 
            topclones_sam2Prod[, indSam2]), ]
        commonClonotypes1Unprod = topclones_sam1Unprod[topclones_sam1Unprod[, 
            indSam1] \%in\% intersect(topclones_sam1Unprod[, indSam1], 
            topclones_sam2Unprod[, indSam2]), ]
        commonClonotypes2Unprod = topclones_sam2Unprod[topclones_sam2Unprod[, 
            indSam2] \%in\% intersect(topclones_sam1Unprod[, indSam1], 
            topclones_sam2Unprod[, indSam2]), ]
        NinsertionsP1 = mean(commonClonotypes1$N2INSERTION + 
            commonClonotypes1$N1INSERTION)
        NinsertionsP2 = mean(commonClonotypes2$N2INSERTION + 
            commonClonotypes2$N1INSERTION)
        NinsertionsUP1 = mean(commonClonotypes1Unprod$N2INSERTION + 
            commonClonotypes1Unprod$N1INSERTION)
        NinsertionsUP2 = mean(commonClonotypes2Unprod$N2INSERTION + 
            commonClonotypes2Unprod$N1INSERTION)
        NinsertionProd1 = c(NinsertionProd1, NinsertionsP1)
        NinsertionProd2 = c(NinsertionProd2, NinsertionsP2)
        NinsertionUnprod1 = c(NinsertionUnprod1, NinsertionsUP1)
        NinsertionUnprod2 = c(NinsertionUnprod2, NinsertionsUP2)
        selectedCols <- c("NUCLEOTIDE", "AMINOACID", "COUNT", 
            "VFAMILYNAME", "VGENENAME", "JGENENAME", "CDR3LENGTH", 
            "SEQUENCESTATUS")
        selectedColsIdx <- which(colnames(commonClonotypes1) \%in\% 
            selectedCols)
        resultDir = paste("Pairwise_Overlap_results_", "Productive_Unproductive(", 
            unProdStatus, ")", sep = "")
        dir.create(file.path(resultDir), showWarnings = FALSE)
        commonsProdForWriting = head(merge(commonClonotypes1[, 
            selectedColsIdx], commonClonotypes2[, selectedColsIdx], 
            by = "NUCLEOTIDE"), 20)
        commonsUnProdForWriting = merge(commonClonotypes1Unprod[, 
            selectedColsIdx], commonClonotypes2Unprod[, selectedColsIdx], 
            by = "NUCLEOTIDE")
        write.table(commonsProdForWriting, file = paste(resultDir, 
            "/", samplesCompared, "_Common_productive_sequences_top20.txt", 
            sep = ""), row.names = F, sep = "\t")
        write.table(commonsUnProdForWriting, file = paste(resultDir, 
            "/", samplesCompared, "_Common_Unproductive_sequences_all.txt", 
            sep = ""), row.names = F, sep = "\t")
        if (nrow(toadd) > 0) {
            newrow = data.frame(overlapValuesProd, overlapValuesTotalSam1Prod, 
                overlapValuesTotalSam2Prod, overlapValuesUnprod, 
                overlapValuesTotalSam1Unprod, overlapValuesTotalSam2Unprod, 
                overlapValuesProdJI, overlapValuesUnprodJI, overlapValuesProdJI_abundancebased, 
                overlapValuesUnprodJI_abundancebased, overlapValuesProdJI_abundancebased2, 
                overlapValuesUnprodJI_abundancebased2, NinsertionProd1, 
                NinsertionProd2, NinsertionUnprod1, NinsertionUnprod2)
            colnames(newrow) <- c("OverlapCountProductive", "OverlapCountTotalProductiveSam1", 
                "OverlapCountTotalProductiveSam2", "OverlapCountUnProductive", 
                "OverlapCountTotalUnProductiveSam1", "OverlapCountTotalUnProductiveSam2", 
                "OverlapProd_JI", "OverlapUnProductive_JI", "AbundanceBasedProd_JI", 
                "AbundanceBasedUnProductive_JI", "AbundanceBasedProd_JI2", 
                "AbundanceBasedUnProductive_JI2", "NinsertionProdSample1", 
                "NinsertionProdSample2", "NinsertionUnProductiveSample1", 
                "NinsertionUnProductiveSample2")
            toadd = rbind(toadd, newrow)
        }
        else {
            toadd = data.frame(overlapValuesProd, overlapValuesTotalSam1Prod, 
                overlapValuesTotalSam2Prod, overlapValuesUnprod, 
                overlapValuesTotalSam1Unprod, overlapValuesTotalSam2Unprod, 
                overlapValuesProdJI, overlapValuesUnprodJI, overlapValuesProdJI_abundancebased, 
                overlapValuesUnprodJI_abundancebased, overlapValuesProdJI_abundancebased2, 
                overlapValuesUnprodJI_abundancebased2, NinsertionProd1, 
                NinsertionProd2, NinsertionUnprod1, NinsertionUnprod2)
            colnames(toadd) <- c("OverlapCountProductive", "OverlapCountTotalProductiveSam1", 
                "OverlapCountTotalProductiveSam2", "OverlapCountUnProductive", 
                "OverlapCountTotalUnProductiveSam1", "OverlapCountTotalUnProductiveSam2", 
                "OverlapProd_JI", "OverlapUnProductive_JI", "AbundanceBasedProd_JI", 
                "AbundanceBasedUnProductive_JI", "AbundanceBasedProd_JI2", 
                "AbundanceBasedUnProductive_JI2", "NinsertionProdSample1", 
                "NinsertionProdSample2", "NinsertionUnProductiveSample1", 
                "NinsertionUnProductiveSample2")
        }
    }
    resultTable = data.frame(m, toadd)
    outputFileName = paste(resultDir, "/", "Pairwise_Overlap_Counts_", 
        seqType, "_Productive_Unproductive(", unProdStatus, ").txt", 
        sep = "")
    write.table(resultTable, file = outputFileName, row.names = T, 
        sep = "\t")
    return(resultTable)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
