% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RepAnFns.R
\name{getOutOfFrame}
\alias{getOutOfFrame}
\title{Get out of frame CDR3s from a repertoire data (Repseq sample) in immunoseq format}
\usage{
getOutOfFrame(rep)
}
\arguments{
\item{rep}{name of the repertoire sample (which has already been read in with for example with read.table)}
}
\value{
a subset of the repertoire (rep) containing only out of frame CDR3s is returned
}
\description{
Get out of frame CDR3s from a repertoire data (Repseq sample) in immunoseq format
}
