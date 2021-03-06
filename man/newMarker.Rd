% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marker.R
\name{newMarker}
\alias{newMarker}
\title{Internal marker constructor}
\usage{
newMarker(
  alleleMatrixInt,
  alleles,
  afreq,
  name = NA_character_,
  chrom = NA_character_,
  posMb = NA_real_,
  mutmod = NULL,
  pedmembers,
  sex
)
}
\arguments{
\item{alleleMatrixInt}{An integer matrix.}

\item{alleles}{A character vector.}

\item{afreq}{A numeric vector.}

\item{name}{A character of length 1.}

\item{chrom}{A character of length 1.}

\item{posMb}{A numeric of length 1.}

\item{mutmod}{A mutation model.}

\item{pedmembers}{A character vector.}

\item{sex}{An integer vector.}
}
\value{
A \code{marker} object.
}
\description{
This is the internal constructor of \code{marker} objects. It does not do any
input validation and should only be used in programming scenarios, and only
if you know what you are doing. Most users are recommended to use the regular
constructor \code{\link[=marker]{marker()}}.
}
\details{
See \code{\link[=marker]{marker()}} for more details about the marker attributes.
}
\examples{

newMarker(matrix(c(1L, 0L, 1L, 1L, 0L, 2L), ncol = 2),
          alleles = c("A", "B"), afreq = c(0.1, 0.9), name = "M",
          pedmembers = c("1", "2", "3"), sex = c(1L, 2L, 1L))

}
