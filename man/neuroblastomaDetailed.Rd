\name{neuroblastomaDetailed}
\alias{neuroblastomaDetailed}
\docType{data}
\title{
Detailed annotations of the neuroblastoma data
}
\description{
  An annotation is the number of breakpoints that an expert
  expects of a segmentation model in a certain region,
  after visual inspection of the
  scatterplot of data.
}
\usage{data(neuroblastomaDetailed)}
\format{
  A data frame with 4359 observations on the following 5 variables.
  \describe{
    \item{\code{profile.id}}{a factor with levels corresponding to the
      profile.id column of \code{neuroblastoma$profiles}.}
    \item{\code{chromosome}}{idem for \code{neuroblastoma$chromosome}.}
    \item{\code{min}}{first position of the annotated region in base pairs.}
    \item{\code{max}}{idem for the last position.}
    \item{\code{annotation}}{factor indicating the number of breakpoints
      in this region:
      \code{>0breakpoints} means at least one breakpoint,
      \code{1breakpoint} means exactly 1 breakpoint,
      \code{normal} means exactly 0 breakpoints.}
  }
}
\details{
The neuroblastoma data are a set of 575 DNA copy number profiles of
neuroblastoma tumors, available as
\code{data(neuroblastoma,package="neuroblastoma")}.
That package provides the "original" set
of up to 6 annotated regions per profile as
\code{neuroblastoma$annotations}. There is at most 1 annotation per
chromosome, and 2 types of annotations: \code{breakpoint} means 1 or more
breakpoints and \code{normal} means exactly 0 breakpoints. These data were
made by Gudrun Schleiermacher and Isabelle Janoueix-Lerosey, by
typing 0 or 1 in a spreadsheet after visual inspection of the profiles.

This package provides a different set
of annotations of the same data. We say they are detailed since there is
often more than 1 annotation per chromosome, and there is another type
of annotation: \code{1breakpoint} means there is exactly 1 breakpoint in that
region. These annotations were created by Toby Dylan Hocking and
Valentina Boeva using GUIs which allow drawing regions on the plotted
data.
}
\source{
http://cbio.ensmp.fr/~thocking/neuroblastoma/annotations.csv
}
\keyword{datasets}
