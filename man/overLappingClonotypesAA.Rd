\name{overLappingClonotypesAA}
\alias{overLappingClonotypesAA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
overLappingClonotypesAA(reps, cutoff = 20, cutoffby = "perc")
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
function (reps, cutoff = 20, cutoffby = "perc") 
{
    cutoffby = toupper(cutoffby)
    m <- as.data.frame(expand.grid(reps, reps))
    colnames(m) <- c("Sample1", "Sample2")
    toadd = data.frame()
    for (r in 1:nrow(m)) {
        overlapValuesProd = c()
        overlapValuesUnprod = c()
        overlapValuesProdJI = c()
        overlapValuesUnprodJI = c()
        overlapValuesProdJI_abundancebased = c()
        overlapValuesUnprodJI_abundancebased = c()
        NinsertionProd1 = c()
        NinsertionUnprod1 = c()
        NinsertionProd2 = c()
        NinsertionUnprod2 = c()
        print(r)
        sam1 = m[r, 1]
        sam2 = m[r, 2]
        sam1 <- get(as.character(sam1))
        sam2 <- get(as.character(sam2))
        sam1Prod <- getProductive(sam1)
        sam1Unprod <- getUnproductive(sam1)
        sam2Prod <- getProductive(sam2)
        sam2Unprod <- getUnproductive(sam2)
        if (cutoffby == "PERC") {
            topclones_sam1Prod = getAbundantClonotypesbyperc(sam1Prod, 
                cutoff)
            topclones_sam1Unprod = getAbundantClonotypesbyperc(sam1Unprod, 
                cutoff)
            topclones_sam2Prod = getAbundantClonotypesbyperc(sam2Prod, 
                cutoff)
            topclones_sam2Unprod = getAbundantClonotypesbyperc(sam2Unprod, 
                cutoff)
            commonClonesProd = length(intersect(topclones_sam1Prod$AMINOACID, 
                topclones_sam2Prod$AMINOACID))
            unionProd = length(union(topclones_sam1Prod$AMINOACID, 
                topclones_sam2Prod$AMINOACID))
            commonClonesUnprod = length(intersect(topclones_sam1Unprod$AMINOACID, 
                topclones_sam2Unprod$AMINOACID))
            unionUnprod = length(union(topclones_sam1Unprod$AMINOACID, 
                topclones_sam2Unprod$AMINOACID))
            if (unionUnprod == 0) {
                unionUnprod = 1
            }
            if (unionProd == 0) {
                unionProd = 1
            }
            overlapValuesProd = c(overlapValuesProd, as.numeric(commonClonesProd))
            overlapValuesUnprod = c(overlapValuesUnprod, as.numeric(commonClonesUnprod))
            overlapValuesProdJI = c(overlapValuesProdJI, commonClonesProd/unionProd)
            overlapValuesUnprodJI = c(overlapValuesUnprodJI, 
                commonClonesUnprod/unionUnprod)
            commonClonesProd_U_data = topclones_sam1Prod[topclones_sam1Prod$AMINOACID \%in\% 
                intersect(topclones_sam1Prod$AMINOACID, topclones_sam2Prod$AMINOACID), 
                ]
            commonClonesProd_U = sum(commonClonesProd_U_data$COUNT)/sum(topclones_sam1Prod$COUNT)
            commonClonesProd_V_data = topclones_sam2Prod[topclones_sam2Prod$AMINOACID \%in\% 
                intersect(topclones_sam1Prod$AMINOACID, topclones_sam2Prod$AMINOACID), 
                ]
            commonClonesProd_V = sum(commonClonesProd_V_data$COUNT)/sum(topclones_sam2Prod$COUNT)
            uv = commonClonesProd_U * commonClonesProd_V
            uplusv = commonClonesProd_U + commonClonesProd_V
            AbundanceBasedJIProd = uv/(uplusv - uv)
            overlapValuesProdJI_abundancebased = c(overlapValuesProdJI_abundancebased, 
                as.numeric(AbundanceBasedJIProd))
            commonClonesUnProd_U_data = topclones_sam1Unprod[topclones_sam1Unprod$AMINOACID \%in\% 
                intersect(topclones_sam1Unprod$AMINOACID, topclones_sam2Unprod$AMINOACID), 
                ]
            commonClonesUnProd_U = sum(commonClonesUnProd_U_data$COUNT)/sum(topclones_sam1Unprod$COUNT)
            commonClonesUnProd_V_data = topclones_sam2Unprod[topclones_sam2Unprod$AMINOACID \%in\% 
                intersect(topclones_sam1Unprod$AMINOACID, topclones_sam2Unprod$AMINOACID), 
                ]
            commonClonesUnProd_V = sum(commonClonesUnProd_V_data$COUNT)/sum(topclones_sam2Unprod$COUNT)
            uv_UnProd = commonClonesUnProd_U * commonClonesUnProd_V
            uplusv_UnProd = commonClonesUnProd_U + commonClonesUnProd_V
            AbundanceBasedJIUnProd = uv_UnProd/(uplusv_UnProd - 
                uv_UnProd)
            overlapValuesUnprodJI_abundancebased = c(overlapValuesUnprodJI_abundancebased, 
                as.numeric(AbundanceBasedJIUnProd))
            commonClonotypes1 = topclones_sam1Prod[topclones_sam1Prod$AMINOACID \%in\% 
                intersect(topclones_sam1Prod$AMINOACID, topclones_sam2Prod$AMINOACID), 
                ]
            commonClonotypes2 = topclones_sam2Prod[topclones_sam2Prod$AMINOACID \%in\% 
                intersect(topclones_sam1Prod$AMINOACID, topclones_sam2Prod$AMINOACID), 
                ]
            commonClonotypes1Unprod = topclones_sam1Unprod[topclones_sam1Unprod$AMINOACID \%in\% 
                intersect(topclones_sam1Unprod$AMINOACID, topclones_sam2Unprod$AMINOACID), 
                ]
            commonClonotypes2Unprod = topclones_sam2Unprod[topclones_sam2Unprod$AMINOACID \%in\% 
                intersect(topclones_sam1Unprod$AMINOACID, topclones_sam2Unprod$AMINOACID), 
                ]
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
            commonClonesProd = length(intersect(topclones_sam1Prod$AMINOACID, 
                topclones_sam2Prod$AMINOACID))
            unionProd = length(union(topclones_sam1Prod$AMINOACID, 
                topclones_sam2Prod$AMINOACID))
            commonClonesUnprod = length(intersect(topclones_sam1Unprod$AMINOACID, 
                topclones_sam2Unprod$AMINOACID))
            unionUnprod = length(union(topclones_sam1Unprod$AMINOACID, 
                topclones_sam2Unprod$AMINOACID))
            if (unionUnprod == 0) {
                unionUnprod = 1
            }
            if (unionProd == 0) {
                unionProd = 1
            }
            overlapValuesProd = c(overlapValuesProd, as.numeric(commonClonesProd))
            overlapValuesUnprod = c(overlapValuesUnprod, as.numeric(commonClonesUnprod))
            overlapValuesProdJI = c(overlapValuesProdJI, commonClonesProd/unionProd)
            overlapValuesUnprodJI = c(overlapValuesUnprodJI, 
                commonClonesUnprod/unionUnprod)
            commonClonesProd_U_data = topclones_sam1Prod[topclones_sam1Prod$AMINOACID \%in\% 
                intersect(topclones_sam1Prod$AMINOACID, topclones_sam2Prod$AMINOACID), 
                ]
            commonClonesProd_U = sum(commonClonesProd_U_data$COUNT)/sum(topclones_sam1Prod$COUNT)
            commonClonesProd_V_data = topclones_sam2Prod[topclones_sam2Prod$AMINOACID \%in\% 
                intersect(topclones_sam1Prod$AMINOACID, topclones_sam2Prod$AMINOACID), 
                ]
            commonClonesProd_V = sum(commonClonesProd_V_data$COUNT)/sum(topclones_sam2Prod$COUNT)
            uv = commonClonesProd_U * commonClonesProd_V
            uplusv = commonClonesProd_U + commonClonesProd_V
            AbundanceBasedJIProd = uv/(uplusv - uv)
            overlapValuesProdJI_abundancebased = c(overlapValuesProdJI_abundancebased, 
                as.numeric(AbundanceBasedJIProd))
            commonClonesUnProd_U_data = topclones_sam1Unprod[topclones_sam1Unprod$AMINOACID \%in\% 
                intersect(topclones_sam1Unprod$AMINOACID, topclones_sam2Unprod$AMINOACID), 
                ]
            commonClonesUnProd_U = sum(commonClonesUnProd_U_data$COUNT)/sum(topclones_sam1Unprod$COUNT)
            commonClonesUnProd_V_data = topclones_sam2Unprod[topclones_sam2Unprod$AMINOACID \%in\% 
                intersect(topclones_sam1Unprod$AMINOACID, topclones_sam2Unprod$AMINOACID), 
                ]
            commonClonesUnProd_V = sum(commonClonesUnProd_V_data$COUNT)/sum(topclones_sam2Unprod$COUNT)
            uv_UnProd = commonClonesUnProd_U * commonClonesUnProd_V
            uplusv_UnProd = commonClonesUnProd_U + commonClonesUnProd_V
            AbundanceBasedJIUnProd = uv_UnProd/(uplusv_UnProd - 
                uv_UnProd)
            overlapValuesUnprodJI_abundancebased = c(overlapValuesUnprodJI_abundancebased, 
                as.numeric(AbundanceBasedJIUnProd))
            commonClonotypes1 = topclones_sam1Prod[topclones_sam1Prod$AMINOACID \%in\% 
                intersect(topclones_sam1Prod$AMINOACID, topclones_sam2Prod$AMINOACID), 
                ]
            commonClonotypes2 = topclones_sam2Prod[topclones_sam2Prod$AMINOACID \%in\% 
                intersect(topclones_sam1Prod$AMINOACID, topclones_sam2Prod$AMINOACID), 
                ]
            commonClonotypes1Unprod = topclones_sam1Unprod[topclones_sam1Unprod$AMINOACID \%in\% 
                intersect(topclones_sam1Unprod$AMINOACID, topclones_sam2Unprod$AMINOACID), 
                ]
            commonClonotypes2Unprod = topclones_sam2Unprod[topclones_sam2Unprod$AMINOACID \%in\% 
                intersect(topclones_sam1Unprod$AMINOACID, topclones_sam2Unprod$AMINOACID), 
                ]
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
        }
        if (nrow(toadd) > 0) {
            newrow = data.frame(overlapValuesProd, overlapValuesUnprod, 
                overlapValuesProdJI, overlapValuesUnprodJI, overlapValuesProdJI_abundancebased, 
                overlapValuesUnprodJI_abundancebased, NinsertionProd1, 
                NinsertionProd2, NinsertionUnprod1, NinsertionUnprod2)
            colnames(newrow) <- c("OverlapCountProductive", "OverlapCountUnProductive", 
                "OverlapProd_JI", "OverlapUnProd_JI", "AbundanceBased_JI", 
                "AbundanceBased_JI", "NinsertionProdSample1", 
                "NinsertionProdSample2", "NinsertionUnProdSample1", 
                "NinsertionUnProdSample2")
            toadd = rbind(toadd, newrow)
        }
        else {
            toadd = data.frame(overlapValuesProd, overlapValuesUnprod, 
                overlapValuesProdJI, overlapValuesUnprodJI, overlapValuesProdJI_abundancebased, 
                overlapValuesUnprodJI_abundancebased, NinsertionProd1, 
                NinsertionProd2, NinsertionUnprod1, NinsertionUnprod2)
            colnames(toadd) <- c("OverlapCountProductive", "OverlapCountUnProductive", 
                "OverlapProd_JI", "OverlapUnProd_JI", "AbundanceBased_JI", 
                "AbundanceBased_JI", "NinsertionProdSample1", 
                "NinsertionProdSample2", "NinsertionUnProdSample1", 
                "NinsertionUnProdSample2")
        }
    }
    print(nrow(toadd))
    resultTable = data.frame(m, toadd)
    forProdselected = resultTable[, c(1, 2, 5)]
    forProdselected = forProdselected[order(forProdselected$Sample1), 
        ]
    l = length(unique(forProdselected$Sample1))
    mat = matrix(as.numeric(forProdselected[, 3]), nrow = l, 
        ncol = l, byrow = T)
    mat[is.nan(mat)] = 0
    colnames(mat) <- unique(forProdselected$Sample1)
    rownames(mat) <- unique(forProdselected$Sample1)
    mat = signif(mat, digits = 2)
    hc = hclust(as.dist(log10(1 - mat)))
    flname = paste("ProdJI_dendrogram.pdf")
    pdf(file = flname, useDingbats = FALSE)
    plot(hc)
    dev.off()
    forUnprodselected = resultTable[, c(1, 2, 6)]
    forUnprodselected = forUnprodselected[order(forUnprodselected$Sample1), 
        ]
    l = length(unique(forUnprodselected$Sample1))
    mat = matrix(as.numeric(forUnprodselected[, 3]), nrow = l, 
        ncol = l, byrow = T)
    mat[is.nan(mat)] = 0
    colnames(mat) <- unique(forUnprodselected$Sample1)
    rownames(mat) <- unique(forUnprodselected$Sample1)
    mat = signif(mat, digits = 2)
    hc = hclust(as.dist(log10(1 - mat)))
    flname = paste("UnprodJI_dendrogram.pdf")
    pdf(file = flname, useDingbats = FALSE)
    plot(hc)
    dev.off()
    write.table(resultTable, file = "Pairwise AA Overlap Counts.txt", 
        row.names = T, sep = "\t")
    return(resultTable)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
