\name{smoothers}
\alias{smoothers}
\docType{data}
\title{smoothers}
\description{This is a list of functions, each of which must return a matrix of
smoothed profiles. The first argument of each function is a
data.frame that represents a copy number profile, with at least
columns: position logratio chromosome. We assume that positions
are already sorted in ascending order p_1 < p_2. The second
argument to each of these functions should be a vector of
smoothing parameters, and there should be a default value. The
matrix returned has 1 row for each parameter, and 1 column for
each position.}
\usage{"smoothers"}





