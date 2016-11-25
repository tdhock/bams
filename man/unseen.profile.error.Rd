\name{unseen.profile.error}
\alias{unseen.profile.error}
\title{unseen profile error}
\description{Do n/t-fold cross-validation to estimate the error of global
models with a small training set.}
\usage{unseen.profile.error(stats, prof.per.train)}
\arguments{
  \item{stats}{list with arrays for errors, false.positive, and false.negative}
  \item{prof.per.train}{t = approximate number of annotated profiles per training set.}
}

\value{matrix of estimated test errors, nfolds x 3}

\author{Toby Dylan Hocking}





