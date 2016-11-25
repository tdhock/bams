\name{exportData}
\alias{exportData}
\title{exportData}
\description{Save profiles and annotations as csv.}
\usage{exportData(data.list = NULL, out.dir = Sys.getenv("HOME"), quote = FALSE, 
    row.names = FALSE, col.names = TRUE, sep = ",", ...)}
\arguments{
  \item{data.list}{List of profiles and annotations to save to disk as csv.}
  \item{out.dir}{Directory where to save them.}
  \item{quote}{passed to write.table}
  \item{row.names}{passed to write.table}
  \item{col.names}{passed to write.table}
  \item{sep}{passed to write.table}
  \item{\dots}{additional arguments for write.table}
}



\author{Toby Dylan Hocking}





