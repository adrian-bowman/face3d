\name{resample.face3d}

\alias{resample.face3d}

\title{Resampling a set of points such that each point remains attached to the surface. }

\description{
  This function resamples a set of points such that each point remains attached to the surface.
}

\usage{
resample.face3d(path, n = 50, method = "spline", threshold = 1e-5)
}

\arguments{
\item{path}{a matrix with three columns containing the co-ordinates of successive points along the path.}
\item{n}{the number of points in the resampled set.}
\item{method}{a character string identifying the method of interpretation, either \code{"linear"} or \code{"spline"}.}
\item{threshold}{the iteration level for the \code{"spline"} method.}
}

\details{
A matrix of points of the specified \code{n} number is returned.
}

\value{
A three-dimensiona matrix of resampled points on the path.
}

