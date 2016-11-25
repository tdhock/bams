\name{run.pelt}
\alias{run.pelt}
\title{run pelt}
\description{Smooth a profile using the PELT algorithm.}
\usage{run.pelt(profile, penalty = "SIC", values = 0, FUN = cpt.mean)}
\arguments{
  \item{profile}{A profile data.frame.}
  \item{penalty}{character specifying the penalty to use.}
  \item{values}{vector of penalty parameters to try.}
  \item{FUN}{PELT function to use.}
}

\value{Matrix of smoothed profiles: nparam x nprobes.}

\author{Toby Dylan Hocking}





