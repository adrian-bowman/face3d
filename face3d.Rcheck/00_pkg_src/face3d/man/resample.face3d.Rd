\name{resample.face3d}

\alias{resample.face3d}

\title{Resampling a set of points such that each point remains attached to the surface. }

\description{
  This function resamples a set of points such that each point remains attached to the surface.
}

\usage{
resample.face3d(shape, n = 50, threshold = 1e-5)
}

\arguments{
\item{shape}{The original matrix of points.}
\item{n}{The chosen number of points in the resampled set.}
\item{threshold}{The iteration level.}

}

\details{
A matrix of points of the specified \code{n} number is returned}

\value{
\item{shaper}{The matrix of points resampled.}
}

