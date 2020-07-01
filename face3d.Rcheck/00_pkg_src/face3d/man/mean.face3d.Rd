\name{mean.face3d}

\alias{mean.face3d}

\title{Finds the mean of a surface}

\description{
  This function finds the mean of a three dimensional surface given associated curvature values.
}

\usage{
mean.face3d(shape, values)}

\arguments{
\item{shape}{Either a \code{face3d object} or a subset thereof.}
\item{values}{Curvature values associated with each point.} 
}

\details{
An overall mean is calculated for the object.
}

\value{
The mean value of the surface is returned.
}

\examples{
\dontrun{
   data(face)
   mean   <- mean.face3d(shape, values)
}}
