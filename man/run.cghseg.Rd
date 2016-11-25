\name{run.cghseg}
\alias{run.cghseg}
\title{run cghseg}
\description{Run cghseg maximum likelihood DP segmentation and pick the model
using picker.}
\usage{run.cghseg(profile, picker)}
\arguments{
  \item{profile}{Profile data.frame.}
  \item{picker}{Function that gets arguments const.lines (a data.frame with one
line for each segment in the maximum likelihood segmentation for a
chrom), smoothed.mat (matrix of smoothed profile for a chrom: kmax
x nprobes), Y (vector of logratio measurements for a chrom), kmax
(maximum number of segments to consider), n (number of probes),
and should return the chosen smoothing vector for the chrom.}
}

\value{Smoothing matrix nparam x nprobes.}

\author{Toby Dylan Hocking}





