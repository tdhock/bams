\name{each.chrom}
\alias{each.chrom}
\title{each chrom}
\description{Apply a smoothing function independently to each chromosome of a
profile.}
\usage{each.chrom(profile, FUN)}
\arguments{
  \item{profile}{Profile data.frame.}
  \item{FUN}{Function that will take a profile data.frame for one chromosome
and return a smoothing matrix for that chromosome: nparam x
nprobes.}
}

\value{Matrix of smoothed profiles for the entire profile: nparam x
nprobes.}

\author{Toby Dylan Hocking}





