\name{getAggregateCount}
\alias{getAggregateCount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getAggregateCount(reps, countfor = "V", seqStatus = "IN", lev = "A", countsorpro = "prop")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reps}{
%%     ~~Describe \code{reps} here~~
}
  \item{countfor}{
%%     ~~Describe \code{countfor} here~~
}
  \item{seqStatus}{
%%     ~~Describe \code{seqStatus} here~~
}
  \item{lev}{
%%     ~~Describe \code{lev} here~~
}
  \item{countsorpro}{
%%     ~~Describe \code{countsorpro} here~~
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
function (reps, countfor = "V", seqStatus = "IN", lev = "A", 
    countsorpro = "prop") 
{
    countfor = toupper(countfor)
    seqStatus = toupper(seqStatus)
    lev = toupper(lev)
    countsorpro = toupper(countsorpro)
    print(countfor)
    if (countfor == "V") {
        features = getAllVgenes(reps)
        m <- as.data.frame(features)
        colnames(m) <- c("VGENENAME")
        for (rep in reps) {
            print(rep)
            p <- get(rep)
            if (seqStatus == "IN") {
                p <- getProductive(p)
            }
            else if (seqStatus == "OUT") {
                p <- getOutOfFrame(p)
            }
            else if (seqStatus == "STOP") {
                p <- getHasStop(p)
            }
            print(nrow(p))
            if (lev == "N") {
                dAgg <- aggregate(NUCLEOTIDE ~ VGENENAME, data = p, 
                  FUN = length)
            }
            else {
                dAgg <- aggregate(AMINOACID ~ VGENENAME, data = p, 
                  FUN = length)
            }
            if (countsorpro == "PROP") {
                dAgg[, 2] <- dAgg[, 2]/sum(dAgg[, 2])
            }
            colnames(dAgg) <- c("VGENENAME", rep)
            m <- merge(m, dAgg, by = "VGENENAME", all.x = T)
        }
        m[is.na(m)] <- 0
        rownames(m) <- m[, 1]
        m <- m[, -1]
    }
    else if (countfor == "J") {
        features = getAllJgenes(reps)
        m <- as.data.frame(features)
        colnames(m) <- c("JGENENAME")
        for (rep in reps) {
            print(rep)
            p <- get(rep)
            if (seqStatus == "IN") {
                p <- getProductive(p)
            }
            else if (seqStatus == "OUT") {
                p <- getOutOfFrame(p)
            }
            else if (seqStatus == "STOP") {
                p <- getHasStop(p)
            }
            print(nrow(p))
            if (lev == "N") {
                dAgg <- aggregate(NUCLEOTIDE ~ JGENENAME, data = p, 
                  FUN = length)
            }
            else {
                dAgg <- aggregate(AMINOACID ~ JGENENAME, data = p, 
                  FUN = length)
            }
            if (countsorpro == "PROP") {
                dAgg[, 2] <- dAgg[, 2]/sum(dAgg[, 2])
            }
            colnames(dAgg) <- c("JGENENAME", rep)
            m <- merge(m, dAgg, by = "JGENENAME", all.x = T)
        }
        m[is.na(m)] <- 0
        rownames(m) <- m[, 1]
        m <- m[, -1]
    }
    else if (countfor == "VF") {
        features = getAllVfamilies(reps)
        m <- as.data.frame(features)
        colnames(m) <- c("VFAMILYNAME")
        for (rep in reps) {
            print(rep)
            p <- get(rep)
            if (seqStatus == "IN") {
                p <- getProductive(p)
            }
            else if (seqStatus == "OUT") {
                p <- getOutOfFrame(p)
            }
            else if (seqStatus == "STOP") {
                p <- getHasStop(p)
            }
            print(nrow(p))
            if (lev == "N") {
                dAgg <- aggregate(NUCLEOTIDE ~ VFAMILYNAME, data = p, 
                  FUN = length)
            }
            else {
                dAgg <- aggregate(AMINOACID ~ VFAMILYNAME, data = p, 
                  FUN = length)
            }
            if (countsorpro == "PROP") {
                dAgg[, 2] <- dAgg[, 2]/sum(dAgg[, 2])
            }
            colnames(dAgg) <- c("VFAMILYNAME", rep)
            m <- merge(m, dAgg, by = "VFAMILYNAME", all.x = T)
        }
        m[is.na(m)] <- 0
        rownames(m) <- m[, 1]
        m <- m[, -1]
    }
    else if (countfor == "VJ") {
        features = getAllVJgenes(reps)
        m <- features
        colnames(m) <- c("VGENENAME", "JGENENAME")
        for (rep in reps) {
            print(rep)
            p <- get(rep)
            if (seqStatus == "IN") {
                p <- getProductive(p)
            }
            else if (seqStatus == "OUT") {
                p <- getOutOfFrame(p)
            }
            else if (seqStatus == "STOP") {
                p <- getHasStop(p)
            }
            print(nrow(p))
            if (lev == "N") {
                dAgg <- aggregate(NUCLEOTIDE ~ VGENENAME + JGENENAME, 
                  data = p, FUN = length)
            }
            else {
                dAgg <- aggregate(AMINOACID ~ VGENENAME + JGENENAME, 
                  data = p, FUN = length)
            }
            if (countsorpro == "PROP") {
                dAgg[, 3] <- dAgg[, 3]/sum(dAgg[, 3])
            }
            colnames(dAgg) <- c("VGENENAME", "JGENENAME", rep)
            m <- merge(m, dAgg, by = c("VGENENAME", "JGENENAME"), 
                all.x = T)
        }
        m[is.na(m)] <- 0
        m[, 1] <- paste(m[, 1], m[, 2])
        m <- m[, -2]
        rownames(m) <- m[, 1]
        m <- m[, -1]
    }
    return(m)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
