\name{run.smoothers}
\alias{run.smoothers}
\title{run smoothers}
\description{Run several smoothers on a profile, quantifying agreement of the
smoother with breakpoint annotations.}
\usage{run.smoothers(profile, breakpoint.labels, smooth.funs = smoothers, 
    loss = several.breakpoints, tosave = c("errors", "seconds", 
        "parameters", "breakpoint.labels"), db = file.path(Sys.getenv("HOME"), 
        "smooth"))}
\arguments{
  \item{profile}{Profile data.frame.}
  \item{breakpoint.labels}{Annotation data.frame.}
  \item{smooth.funs}{List of smoothing functions to apply to the profile.}
  \item{loss}{Loss function of (breakpoint counts,labels) --- both matrices
nparam x nlabels, and should return numeric matrix of same size.}
  \item{tosave}{Variables to save to the db directory.}
  \item{db}{Location to save gzipped result files.}
}

\value{Nothing, the results are saved to files.}

\author{Toby Dylan Hocking}





