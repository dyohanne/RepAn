\name{getClusterMatches}
\alias{getClusterMatches}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getClusterMatches(repSeqObj, matchingMethod = c("hc", "km", "og"), distMethod = "euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repSeqObj}{
%%     ~~Describe \code{repSeqObj} here~~
}
  \item{matchingMethod}{
%%     ~~Describe \code{matchingMethod} here~~
}
  \item{distMethod}{
%%     ~~Describe \code{distMethod} here~~
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
function (repSeqObj, matchingMethod = c("hc", "km", "og"), distMethod = "euclidean") 
{
    if (matchingMethod == "og") {
        matchTable = NULL
        clMatches = c()
        samClusterCentroids = repSeqObj$withinSampleClusters[[1]]$clusterCenters
        subRepertoireNames = 1:nrow(samClusterCentroids)
        for (i in 1:nrow(samClusterCentroids)) {
            for (s2 in 2:length(repSeqObj$samNames)) {
                combinedsams1 <- rbind(samClusterCentroids[i, 
                  ], repSeqObj$withinSampleClusters[[s2]]$clusterCenters)
                imin <- getMatchingCluster(combinedsams1)
                if (!is.null(matchTable) & (imin \%in\% matchTable[, 
                  s2])) {
                  cimins <- c()
                  while (imin \%in\% matchTable[, s2]) {
                    cimins <- c(cimins, imin)
                    combinedsams1 <- rbind(samClusterCentroids[i, 
                      ], repSeqObj$withinSampleClusters[[s2]]$clusterCenters[!rownames(repSeqObj$withinSampleClusters[[s2]]$clusterCenters) \%in\% 
                      as.numeric(cimins), ])
                    if ((length(subRepertoireNames) - length(cimins)) > 
                      1) {
                      imin <- getMatchingCluster(combinedsams1)
                    }
                    else {
                      imin <- subRepertoireNames[!(subRepertoireNames \%in\% 
                        cimins)]
                    }
                  }
                }
                clMatches <- c(clMatches, as.numeric(imin))
            }
            matchTable <- rbind(matchTable, c(i, clMatches))
            clMatches = c()
        }
        colnames(matchTable) <- repSeqObj$samNames
    }
    else {
        combinedCentroids = NULL
        for (sam in repSeqObj$samNames) {
            samClusterCentroids <- repSeqObj$withinSampleClusters[[sam]]$clusterCenters
            ntimes <- nrow(samClusterCentroids)
            rownames(samClusterCentroids) <- paste(rep(sam, ntimes), 
                "_", rownames(samClusterCentroids), sep = "")
            combinedCentroids <- rbind(combinedCentroids, samClusterCentroids)
        }
        if (matchingMethod == "km") {
            maxKacrossSamples = max(sapply(repSeqObj$samNames, 
                function(x) length(table(repSeqObj$withinSampleClusters[[x]]$clusters))))
            minKacrossSamples = min(sapply(repSeqObj$samNames, 
                function(x) length(table(repSeqObj$withinSampleClusters[[x]]$clusters))))
            sils <- c()
            combinedCentroidsdist <- dist(combinedCentroids)
            for (j in minKacrossSamples:maxKacrossSamples) {
                centroidKmcls <- kmeans(combinedCentroids, maxKacrossSamples, 
                  iter.max = 50, nstart = 50)
                sicentroid <- silhouette(centroidKmcls$cluster, 
                  combinedCentroidsdist)
                avgSilWidth <- summary(sicentroid)$avg.width
                nClusSilWidthAboveAverage <- sum(summary(sicentroid)$clus.avg.widths > 
                  avgSilWidth)/j
                compressionGain <- j/nrow(combinedCentroids)
                avgSil <- (avgSilWidth + 1)/2
                kscore <- (nClusSilWidthAboveAverage + avgSil)/2
                sils <- c(sils, avgSilWidth)
            }
            goodK = minKacrossSamples + which.max(sils) - 1
            centroidKmcls <- kmeans(combinedCentroids, goodK, 
                iter.max = 100, nstart = 50)
        }
        else if (matchingMethod == "hc") {
            centroidKmcls <- list()
            dcalc1 = dist.matrix(combinedCentroids, method = distMethod, 
                convert = T, as.dist = TRUE)
            hc <- hclust(dcalc1, method = "complete")
            cls <- cutreeDynamic(hc, distM = as.matrix(dcalc1), 
                minClusterSize = 1, verbose = 0)
            centroidKmcls$cluster <- cls
            names(centroidKmcls$cluster) <- rownames(combinedCentroids)
            goodK <- length(unique(centroidKmcls$cluster))
        }
        matchTable = matrix(data = NA, nrow = goodK, ncol = length(repSeqObj$samNames), 
            dimnames = list(1:goodK, repSeqObj$samNames))
        for (i in 1:length(table(centroidKmcls$cluster))) {
            tempCentroidClustMem <- centroidKmcls$cluster[centroidKmcls$cluster == 
                i]
            availSamples <- sapply(names(tempCentroidClustMem), 
                function(x) strsplit(x, "_")[[1]][1])
            matchingClustersInSamples <- sapply(names(tempCentroidClustMem), 
                function(x) strsplit(x, "_")[[1]][2])
            removedLables <- c()
            if (sum(table(availSamples) > 1) > 0) {
                for (name in names(which(table(availSamples) > 
                  1))) {
                  clustersToBeMerged <- matchingClustersInSamples[availSamples == 
                    name]
                  idxTobeRemoved <- which(availSamples == name)[-1]
                  repSeqObj$withinSampleClusters[[name]]$clusters[repSeqObj$withinSampleClusters[[name]]$clusters \%in\% 
                    clustersToBeMerged] <- as.character(clustersToBeMerged[1])
                  updatedCentroid <- getCenters(repSeqObj$withinSampleClusters[[name]]$seqmers[repSeqObj$withinSampleClusters[[name]]$clusters \%in\% 
                    clustersToBeMerged[1], ], repSeqObj$withinSampleClusters[[name]]$clusters[repSeqObj$withinSampleClusters[[name]]$clusters \%in\% 
                    clustersToBeMerged[1]])
                  repSeqObj$withinSampleClusters[[name]]$clusterCenters[as.character(clustersToBeMerged[1]), 
                    ] <- updatedCentroid
                  removedLables <- c(removedLables, idxTobeRemoved)
                }
            }
            if (length(removedLables) > 0) {
                matchingClustersInSamples <- matchingClustersInSamples[-removedLables]
                availSamples <- availSamples[-removedLables]
            }
            matchTable[i, match(availSamples, repSeqObj$samNames)] <- as.numeric(matchingClustersInSamples)
        }
    }
    return(matchTable)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
