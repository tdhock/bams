\name{estimate.test.error}
\alias{estimate.test.error}
\title{estimate test error}
\description{Do leave-one-out cross-validation on chromosome arms.}
\usage{estimate.test.error(stats)}
\arguments{
  \item{stats}{Named list with arrays errors, false.positive, false.negative,
each of dim nparam x nprof x nfolds.}
}

\value{Named list with elements local, hybrid, global, each a 3 x nfolds
matrix.}

\author{Toby Dylan Hocking}





