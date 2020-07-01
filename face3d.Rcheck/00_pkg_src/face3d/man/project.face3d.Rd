\name{project.face3d}

\alias{project.face3d}

\title{Project points onto the face using supplied normal vectors.}

\description{
This function projects the points \code{x} onto a \code{face3d} object, by moving along the directions specified in the corresponding \code{normals} vectors.
}

\usage{
project.face3d(x, normals, shape)
}

\arguments{
\item{shape}{a \code{face3d} object.}
}

\details{
This function projects the points \code{x} onto a \code{face3d} object, by moving along the directions specified in the corresponding \code{normals} vectors.
}

\value{
A matrix whose rows give the locations of the project points.
}

\examples{
\dontrun{
   # To be written.
}}
